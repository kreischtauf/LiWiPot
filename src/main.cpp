/***************************************************************************
                               main.cpp
                         -------------------
    begin                : Mon Mar 17 2014
    copyright            : (C) 2014 by Christof Kraus
    email                : christof.kraus-csrst@my.mail.de
***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
#include <iostream>
#include <sstream>
#include <fstream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_errno.h>
#include <boost/regex.hpp>

#include "tclap/CmdLine.h"

#include "defs.hh"
#include "PlotBinaryVTK.hh"
#include "BinaryVTKFile.hh"
#include "FieldPatch.hh"
#include "Communicator.hh"
#include "Physics.hh"

#define LIMIT_INTEGRATION_SIZE 100
#define N_INTEG_SAMPLES GSL_INTEG_GAUSS15
std::ofstream dbg;

#define DBGOUT dbg << __LINE__ << ":\t"

double intExBendAcc(double x, void* params);
double intEyBendAcc(double y, void* params);
double intExBendVel(double x, void* params);
double intEyBendVel(double y, void* params);

double ExBendAcc(double y, void* params);
double EyBendAcc(double x, void* params);
double ExBendVel(double y, void* params);
double EyBendVel(double x, void* params);

void errorHandler (const char * reason,
                   const char * file,
                   int line,
                   int gsl_errno);

int main(int argc, char *argv[]){
    Ippl ippl(argc, argv);
    const double lightSpeedInv = 1.0 / Physics::c;
    gsl_set_error_handler(&errorHandler);

    Inform msg(argv[0]);
    Inform msg2all(argv[0], INFORM_ALL_NODES);

    std::ios_base::fmtflags ff = std::cout.flags();
    std::streamsize ss = std::cout.precision();
    IpplTimings::TimerRef mainTimer = IpplTimings::getTimer("Total");
    IpplTimings::TimerRef crunchTimer = IpplTimings::getTimer("crunching");
    IpplTimings::startTimer(mainTimer);

    TCLAP::CmdLine cmdl("lw, an application to calculate the Lienard-Wiechert potentials inside and after a bend", ' ', "0.9");
    TCLAP::ValueArg<unsigned> numSteps  ("t", "numsteps", "number of steps inside the bend", false, 1000, "Integer", cmdl);
    TCLAP::ValueArg<unsigned> numCopies ("N", "numcopies", "number of image charge copies drift", false, 0, "Integer", cmdl);
    TCLAP::ValueArg<unsigned> numCopies2("n", "numcopies2", "number of image charge copies bend", false, 0, "Integer", cmdl);
    TCLAP::ValueArg<double> askedDx     ("x", "dx", "length of cell in x-direction", false, 9.77517e-5, "Double", cmdl);
    TCLAP::ValueArg<double> askedDy     ("y", "dy", "length of cell in y-direction", false, 0.000393701, "Double", cmdl);
    TCLAP::ValueArg<double> phiBend     ("p", "phi","angle of deflection", false, 5.0, "Double", cmdl);
    TCLAP::ValueArg<double> lenAfter    ("l", "length","length of drift after bend", false, 0.10, "Double", cmdl);
    TCLAP::ValueArg<double> width       ("w", "width","width of vacuum pipe", false, 0.10, "Double", cmdl);
    TCLAP::ValueArg<double> ekin        ("E", "Ekin","kinetic energy [MeV]", false, 150.0, "Double", cmdl);
    TCLAP::ValueArg<double> bzext       ("B", "Bz","transverse magnetic field [T]", false, 0.502, "Double", cmdl);
    TCLAP::ValueArg<double> tcharge     ("Q", "Qtotal","total charge [nC]", false, 1.0, "Double", cmdl);
    TCLAP::ValueArg<double> sigmaX      ("s", "sigma_x", "standard deviation of initial distribution", false, 0.0007, "Double", cmdl);
    TCLAP::SwitchArg VTKoutput          ("v", "visual", "produce a VTK output", cmdl, false);
    TCLAP::SwitchArg centeredVTK        ("c", "centered", "lay origin of VTK output at particle position", cmdl, false);
    TCLAP::SwitchArg withoutInitials    ("I", "no_initials", "don't calculate the initial conditions of the em fields", cmdl, false);
    try {
	cmdl.parse(argc, argv);
    } catch(std::exception &e) {
	std::cout << e.what() << std::endl;
	return EXIT_FAILURE;
    }

    const unsigned int Nt = numSteps.getValue();

    const double dx = askedDx.getValue();
    const double dy = askedDy.getValue();
    const double convFactor = std::sqrt(2.0) * sigmaX.getValue(), convFactorInv = 1.0 / convFactor;

    char buffer[18];
    sprintf(buffer, "debug_%05d.txt", Ippl::myNode());
    dbg.open(buffer);

    const double L = lenAfter.getValue();
    const double bz = bzext.getValue();
    const double phi = std::abs(phiBend.getValue()) * bz / std::abs(bz) * Physics::pi / 180.0;
    const double mass = Physics::m_e * 1000;
    const double radius = std::sqrt(ekin.getValue()*ekin.getValue() + 2.0 * ekin.getValue() * mass) * 1e6 / (Physics::c * std::abs(bz));
    const double vwidth = width.getValue();
    const double radiusInv = 1.0 / radius;

    const double gamma = ekin.getValue() / mass + 1.0;
    const double beta = std::sqrt(1.0 - 1.0 / (gamma * gamma));
    const double betasqr = beta*beta;
    const double fourpieps = 0.25  / (Physics::pi * Physics::epsilon_0);

    const Vector_t bunchPosition(radius * sin(phi) + L * cos(phi),
                                 radius * (1.0 - cos(phi)) + L * sin(phi));

    const int Nx = static_cast<int>(2 * std::floor(4 * sigmaX.getValue() / dx + 0.5));
    const int Ny = static_cast<int>(2 * std::floor(4 * sigmaX.getValue() / dy + 0.5));
    const int NyRecv = VTKoutput.getValue()? Ny:static_cast<int>(std::ceil(sin(phi) * dx / dy * Nx) + 1);

    e_dim_tag decomp[] = {PARALLEL, SERIAL, SERIAL};
    double spacing[] = {dx, dy};
    Vector_t origin(bunchPosition(0) - Nx * dx,
                    bunchPosition(1) - NyRecv * dy);
    Vector_t auxOrigin(bunchPosition(0) - (1.5 * Nx - 0.5) * dx,
                       bunchPosition(1) - (NyRecv + 0.5 * Ny - 0.5) * dy);

    Mesh_t mesh(Index(2 * Nx),
                Index(2 * NyRecv),
                spacing,
                origin);
    Mesh_t mesh2(Index(3 * Nx),
                 Index(2 * NyRecv + Ny),
                 spacing,
                 auxOrigin);

    msg << "\n" << std::setprecision(4) << std::fixed
        << "Number of steps:         " << numSteps.getValue()  << "\n"
        << "Bend angle:              " << phiBend.getValue()   << " deg\n"
        << "                         " << phi                  << " rad\n"
        << "Length after bend:       " << lenAfter.getValue()  << " m\n"
        << "Width of pipe:           " << width.getValue()     << " m\n"
        << "Number of image charges: " << numCopies.getValue() << "\n"
        << "Kinetic energy:          " << ekin.getValue()      << " MeV\n"
        << "Magnetic field:          " << bzext.getValue()     << " T\n"
        << "Bend radius:             " << radius               << " m\n"
        << "Total charge:            " << tcharge.getValue()   << " nC\n"
        << "Sigma distribution:      " << sigmaX.getValue()    << " m\n"
        << "\n"
        << "Mesh dimensions:\n";
    std::cout.flags(ff);
    std::cout.precision(ss);

    msg << std::setw(18) << 2 * Nx << " x "
        << std::setw(15) << 2 * NyRecv << "\n"
        << std::setw(18) << std::setprecision(8) << dx << " m "
        << std::setw(15) << std::setprecision(8) << dy << " m\n"
        << "\n"
        << "Source mesh dimensions:\n"
        << std::setw(18) << Nx << " x "
        << std::setw(15) << Ny << "\n"
        << std::setw(18) << std::setprecision(8) << dx << " m "
        << std::setw(15) << std::setprecision(8) << dy << " m\n"
        << "\n"
        << "Computing domain:\n"
        << std::setw(18) << std::setprecision(8) << origin(1) + 2 * NyRecv * dy << "\n"
        << std::setw(18) << std::setprecision(8) << origin(1) << "\n"
        << std::setw(36) << std::setprecision(8) << origin(0)
        << std::setw(18) << std::setprecision(8) << origin(0) + 2 * Nx * dx
        << endl;

    FieldLayout_t FL (mesh, decomp);
    FieldLayout_t FL2(mesh2, decomp);

    FieldPatch<WrappedVector> accAFD(FL2.getLocalNDIndex());
    FieldPatch<WrappedVector> velAFD(FL2.getLocalNDIndex());

    VField_t accEFD(mesh, FL, GCS_t(GUARDCELLSIZE));
    VField_t velEFD(mesh, FL, GCS_t(GUARDCELLSIZE));
    VField_t driftEFD(mesh, FL, GCS_t(GUARDCELLSIZE));
    accEFD = 0.0;
    velEFD = 0.0;
    driftEFD = 0.0;

    std::vector<NDIndex<DIM> > localEFDomains, localAFDomains;
    Utils::getLocalDomains(FL, localEFDomains);
    Utils::getLocalDomains(FL2, localAFDomains);
    for (NDIndex<DIM> & ldom: localEFDomains) {
        ldom[0] = Index(ldom[0].first(), ldom[0].last() + Nx);
        ldom[1] = Index(ldom[1].first(), ldom[1].last() + Ny);
    }

    const double timeToBend = ekin.getValue() > 0.0? radius * phi / (Physics::c * beta): 0.0;
    const double endTime = ekin.getValue() > 0.0? timeToBend + L / (Physics::c * beta): 0.0;
    const double dt = timeToBend / (Nt - 1);

    NDIndex<DIM> lDom = FL.getLocalNDIndex(), lDom2 = FL2.getLocalNDIndex();
    NDIndex<DIM> elem, elem2;

    const double gammasqr = gamma*gamma;
    const Vector3D_t ez(0.0, 0.0, 1.0);

    Vector3D_t X(0.0), V(beta, 0, 0), A(0, betasqr * Physics::c * radiusInv, 0);
    Vector3D_t n;
    Vector3D_t R;

    double charge = 1.0;

    double loadWeightBefore = 1.0;
    double loadWeightBend = 4.0 * 15;
    double loadWeightAfter = 2.0 * 15;
    double totalLoad = ((loadWeightBefore + loadWeightAfter) * (2 * numCopies.getValue() + 1) +
                        loadWeightBend * (numCopies2.getValue() + 1) * (Nt - 1));
    if (withoutInitials.getValue()) totalLoad -= loadWeightBefore * (2 * numCopies.getValue() + 1);
    if (L <= 0.0) totalLoad -= loadWeightAfter * (2 * numCopies.getValue() + 1);

    double totalLoadInv = 1.0 / totalLoad;
    double finishedLoad = 0.0;

    gsl_integration_workspace *w = gsl_integration_workspace_alloc(LIMIT_INTEGRATION_SIZE);
    gsl_integration_workspace *v = gsl_integration_workspace_alloc(LIMIT_INTEGRATION_SIZE);

    IpplTimings::startTimer(crunchTimer);
    if (!withoutInitials.getValue()) {
        const double velPrefactor = fourpieps / (dx * dy);
        const double mincdt = Physics::c * endTime;
        const double maxcdt = 1e17; // just any large number will do

        X = beta * mincdt * Vector3D_t(1.0, 0.0, 0.0);
        V = beta * Vector3D_t(1.0, 0.0, 0.0);
        A = 0.0;
        velAFD = Vector_t(0.0);

        // gsl_function exVel;
        // exVel.function = &intExBendVel;
        gsl_function eyVel;
        eyVel.function = &intEyBendVel;

        std::vector<char> cparams;
        cparams.reserve(10 * sizeof(double) + sizeof(char*));
        double zero = 0.0;
        const char* buffer = reinterpret_cast<const char*>(&zero);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&V(0));
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&V(1));
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&A(0));
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&A(1));
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&gammasqr);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&mincdt);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&maxcdt);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&zero);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&zero);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&v);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(char*));

        double chargeSign = -1.0;
        for (unsigned int nC = 0; nC <= numCopies.getValue(); ++ nC) {
            chargeSign = -chargeSign;
            double chargeDist = nC * vwidth;

            const unsigned int numSides = nC > 0? 2: 1;
            for (unsigned int side = 0; side < numSides; ++ side) {
                if (side == 1) chargeDist = -chargeDist;
                finishedLoad += loadWeightBefore;
                if (std::floor((100.0 * finishedLoad) * totalLoadInv) !=
                    std::floor((100.0 * (finishedLoad - loadWeightBefore)) * totalLoadInv)) {
                    msg << "processing... [" << std::setw(3) << std::setprecision(0) << std::floor((100.0 * finishedLoad) * totalLoadInv) << "%]; " << endl;
                }

                R(2) = 0.0;
                for (int j = lDom2[1].first(); j <= lDom2[1].last(); ++ j) {
                    R(1) = auxOrigin(1) - X(1) + j * dy + chargeDist;

                    for (int i = lDom2[0].first(); i <= lDom2[0].last(); ++ i) {
                        R(0) = auxOrigin(0) - X(0) + i * dx;

                        double intoneXInt = 0.0;
                        double intoneYInt = 0.0;

                        {
                            double lr[] = {R(0) - 0.5 * dx, R(1) - dy};
                            double ur[] = {R(0) + 0.5 * dx, R(1)};
                            *reinterpret_cast<double*>(&cparams[0] + 8 * sizeof(double)) = lr[1];
                            *reinterpret_cast<double*>(&cparams[0] + 9 * sizeof(double)) = ur[1];

                            exVel.params = &cparams[0];

                            double abserr = 0.0;
                            gsl_integration_qag(&exVel,
                                                lr[0],
                                                ur[0],
                                                1e-6,
                                                1e-3,
                                                LIMIT_INTEGRATION_SIZE,
                                                N_INTEG_SAMPLES,
                                                w,
                                                &intoneXInt,
                                                &abserr);
                        }

                        {
                            double lr[] = {R(0) - dx, R(1) - 0.5 * dy};
                            double ur[] = {R(0)     , R(1) + 0.5 * dy};
                            *reinterpret_cast<double*>(&cparams[0] + 8 * sizeof(double)) = lr[0];
                            *reinterpret_cast<double*>(&cparams[0] + 9 * sizeof(double)) = ur[0];

                            eyVel.params = &cparams[0];

                            double abserr = 0.0;
                            gsl_integration_qag(&eyVel,
                                                lr[1],
                                                ur[1],
                                                1e-6,
                                                1e-3,
                                                LIMIT_INTEGRATION_SIZE,
                                                N_INTEG_SAMPLES,
                                                w,
                                                &intoneYInt,
                                                &abserr);
                        }

                        velAFD(i,j) += chargeSign * Vector_t(intoneXInt, intoneYInt);
                    }
                }
            }
        }

        Communicator::communicateFields(velAFD, localAFDomains, localEFDomains, true);

        for (int jo = 0; jo <= Ny; ++ jo) {
            Vector_t Pos(0.0, (jo - Ny / 2) * dy);
            double chargeY = -0.5 * (std::erf((Pos(1) + 0.5 * dy) * convFactorInv) -
                                     std::erf((Pos(1) - 0.5 * dy) * convFactorInv)) * tcharge.getValue() * 1e-9;

            for (int io = 0; io <= Nx; ++ io) {
                Pos(0) = (io - Nx / 2) * dx;

                charge = (0.5 * chargeSign *
                          (std::erf((Pos(0) + 0.5 * dx) * convFactorInv) -
                           std::erf((Pos(0) - 0.5 * dx) * convFactorInv)) * chargeY);

                for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
                    elem[1] = Index(j,j);

                    for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
                        elem[0] = Index(i,i);

                        driftEFD.localElement(elem) += charge * velPrefactor * velAFD(i - io + Nx, j - jo + Ny);
                    }
                }
            }
        }
    }

    if (Nt > 1) {
        const double dphi = phi / (Nt - 1);

        velAFD = Vector_t(0.0);
        accAFD = Vector_t(0.0);

        gsl_function exAcc, eyAcc, exVel, eyVel;
        exAcc.function = &intExBendAcc;
        eyAcc.function = &intEyBendAcc;
        exVel.function = &intExBendVel;
        eyVel.function = &intEyBendVel;

        const double velPrefactor = fourpieps / (dx * dy);
        const double accPrefactor = gammasqr*gammasqr * velPrefactor * lightSpeedInv;

        for (unsigned int k = 0; k < Nt - 1; ++ k) {
            double chargeSign = -1.0;

            const double maxcdt = Physics::c * (endTime - k * dt);
            const double mincdt = maxcdt - Physics::c * dt;

            const double cosphi = cos(k * dphi), sinphi = sin(k * dphi);

            X = radius * Vector3D_t(sinphi, 1.0 - cosphi, 0.0);
            V = beta * Vector3D_t(cosphi, sinphi, 0.0);
            A = beta*beta * Physics::c / radius * Vector3D_t(-sinphi, cosphi, 0.0);

            std::vector<char> cparams;
            cparams.reserve(10 * sizeof(double) + sizeof(char*));
            double zero = 0.0;
            const char* buffer = reinterpret_cast<const char*>(&zero);
            cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
            buffer = reinterpret_cast<const char*>(&V(0));
            cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
            buffer = reinterpret_cast<const char*>(&V(1));
            cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
            buffer = reinterpret_cast<const char*>(&A(0));
            cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
            buffer = reinterpret_cast<const char*>(&A(1));
            cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
            buffer = reinterpret_cast<const char*>(&gammasqr);
            cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
            buffer = reinterpret_cast<const char*>(&mincdt);
            cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
            buffer = reinterpret_cast<const char*>(&maxcdt);
            cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
            buffer = reinterpret_cast<const char*>(&zero);
            cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
            buffer = reinterpret_cast<const char*>(&zero);
            cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
            buffer = reinterpret_cast<const char*>(&v);
            cparams.insert(cparams.end(), buffer, buffer + sizeof(char*));

            for (unsigned int nC = 0; nC <= numCopies2.getValue(); ++ nC) {
                chargeSign = -chargeSign;
                double chargeDist = nC * vwidth;

                finishedLoad += loadWeightBend;
                if (std::floor((100.0 * finishedLoad) * totalLoadInv) !=
                    std::floor((100.0 * (finishedLoad - loadWeightBend)) * totalLoadInv)) {
                    msg << "processing... [" << std::setw(3) << std::setprecision(0)  << std::floor((100.0 * finishedLoad) * totalLoadInv) << "%]; " << endl;
                }

                for (int j = lDom2[1].first(); j <= lDom2[1].last(); ++ j) {
                    R(1) = auxOrigin(1) - X(1) - V(1) * maxcdt + j * dy + cosphi * chargeDist;

                    for (int i = lDom2[0].first(); i <= lDom2[0].last(); ++ i) {
                        R(0) = auxOrigin(0) - X(0) - V(0) * maxcdt + i * dx - sinphi * chargeDist;

                        double exAccInt = 0.0;
                        double eyAccInt = 0.0;
                        double intoneXInt = 0.0;
                        double intoneYInt = 0.0;

                        {
                            double lr[] = {R(0) - 0.5 * dx, R(1) - dy};
                            double ur[] = {R(0) + 0.5 * dx, R(1)};
                            *reinterpret_cast<double*>(&cparams[0] + 8 * sizeof(double)) = lr[1];
                            *reinterpret_cast<double*>(&cparams[0] + 9 * sizeof(double)) = ur[1];

                            exAcc.params = &cparams[0];
                            exVel.params = &cparams[0];

                            double abserr = 0.0;
                            gsl_integration_qag(&exAcc,
                                                lr[0],
                                                ur[0],
                                                1e-6,
                                                1e-3,
                                                LIMIT_INTEGRATION_SIZE,
                                                N_INTEG_SAMPLES,
                                                w,
                                                &exAccInt,
                                                &abserr);

                            gsl_integration_qag(&exVel,
                                                lr[0],
                                                ur[0],
                                                1e-6,
                                                1e-3,
                                                LIMIT_INTEGRATION_SIZE,
                                                N_INTEG_SAMPLES,
                                                w,
                                                &intoneXInt,
                                                &abserr);
                        }

                        {
                            double lr[] = {R(0) - dx, R(1) - 0.5 * dy};
                            double ur[] = {R(0)     , R(1) + 0.5 * dy};
                            *reinterpret_cast<double*>(&cparams[0] + 8 * sizeof(double)) = lr[0];
                            *reinterpret_cast<double*>(&cparams[0] + 9 * sizeof(double)) = ur[0];

                            eyAcc.params = &cparams[0];
                            eyVel.params = &cparams[0];

                            double abserr = 0.0;
                            gsl_integration_qag(&eyAcc,
                                                lr[1],
                                                ur[1],
                                                1e-6,
                                                1e-3,
                                                LIMIT_INTEGRATION_SIZE,
                                                N_INTEG_SAMPLES,
                                                w,
                                                &eyAccInt,
                                                &abserr);

                            gsl_integration_qag(&eyVel,
                                                lr[1],
                                                ur[1],
                                                1e-6,
                                                1e-3,
                                                LIMIT_INTEGRATION_SIZE,
                                                N_INTEG_SAMPLES,
                                                w,
                                                &intoneYInt,
                                                &abserr);
                        }

                        velAFD(i,j) += chargeSign * Vector_t(intoneXInt, intoneYInt);
                        accAFD(i,j) += chargeSign * Vector_t(exAccInt, eyAccInt);
                    }
                }
            }
        }

        Communicator::communicateFields(velAFD, localAFDomains, localEFDomains, true);
        Communicator::communicateFields(accAFD, localAFDomains, localEFDomains, true);

        for (int jo = 0; jo <= Ny; ++ jo) {
            Vector_t Pos(0.0, (jo - Ny / 2) * dy);
            double chargeY = -0.5 * (std::erf((Pos(1) + 0.5 * dy) * convFactorInv) -
                                     std::erf((Pos(1) - 0.5 * dy) * convFactorInv)) * tcharge.getValue() * 1e-9;

            for (int io = 0; io <= Nx; ++ io) {
                Pos(0) = (io - Nx / 2) * dx;

                charge = 0.5 * (std::erf((Pos(0) + 0.5 * dx) * convFactorInv) -
                                std::erf((Pos(0) - 0.5 * dx) * convFactorInv)) * chargeY;


                for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
                    elem[1] = Index(j,j);

                    for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
                        elem[0] = Index(i,i);

                        velEFD.localElement(elem) += charge * velPrefactor * velAFD(i - io + Nx, j - jo + Ny);
                        accEFD.localElement(elem) += charge * accPrefactor * accAFD(i - io + Nx, j - jo + Ny);
                    }
                }
            }
        }
    }

    if (L > 0.0) {
        const double velPrefactor = fourpieps / (dx * dy);
        const double cosphi = cos(phi), sinphi = sin(phi);
        const double maxcdt = Physics::c * (endTime - timeToBend);
        const double mincdt = Physics::c * dt;

        V = beta * Vector3D_t(cosphi, sinphi, 0.0);
        X = Vector3D_t(bunchPosition(0), bunchPosition(1), 0.0);

        velAFD = Vector_t(0.0);

        std::vector<char> cparams;
        cparams.reserve(10 * sizeof(double) + sizeof(char*));
        double zero = 0.0;
        const char* buffer = reinterpret_cast<const char*>(&zero);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&V(0));
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&V(1));
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&zero);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&zero);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&gammasqr);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&mincdt);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&maxcdt);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&zero);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&zero);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(double));
        buffer = reinterpret_cast<const char*>(&v);
        cparams.insert(cparams.end(), buffer, buffer + sizeof(char*));

        gsl_function exVel, eyVel;
        exVel.function = &intExBendVel;
        eyVel.function = &intEyBendVel;

        double chargeSign = -1.0;
        for (unsigned int nC = 0; nC <= numCopies.getValue(); ++ nC) {
            chargeSign = -chargeSign;
            double chargeDist = nC * vwidth;

            const unsigned int numSides = nC > 0? 2: 1;
            for (unsigned int side = 0; side < numSides; ++ side) {
                if (side == 1) chargeDist = -chargeDist;
                finishedLoad += loadWeightAfter;
                if (std::floor((100.0 * finishedLoad) * totalLoadInv) !=
                    std::floor((100.0 * (finishedLoad - loadWeightAfter)) * totalLoadInv)) {
                    msg << "processing... [" << std::setw(3) << std::setprecision(0)
                        << (100.0 * finishedLoad) * totalLoadInv << "%]; " << endl;
                }

                for (int j = lDom2[1].first(); j <= lDom2[1].last(); ++ j) {
                    R(1) = auxOrigin(1) - X(1) + j * dy + chargeDist;

                    for (int i = lDom2[0].first(); i <= lDom2[0].last(); ++ i) {
                        R(0) = auxOrigin(0) - X(0) + i * dx;

                        double intoneXInt = 0.0;
                        double intoneYInt = 0.0;

                        {
                            double lr[] = {R(0) - 0.5 * dx, R(1) - dy};
                            double ur[] = {R(0) + 0.5 * dx, R(1)};
                            *reinterpret_cast<double*>(&cparams[0] + 8 * sizeof(double)) = lr[1];
                            *reinterpret_cast<double*>(&cparams[0] + 9 * sizeof(double)) = ur[1];

                            exVel.params = &cparams[0];

                            double abserr = 0.0;
                            gsl_integration_qag(&exVel,
                                                lr[0],
                                                ur[0],
                                                1e-6,
                                                1e-3,
                                                LIMIT_INTEGRATION_SIZE,
                                                N_INTEG_SAMPLES,
                                                w,
                                                &intoneXInt,
                                                &abserr);
                        }

                        {
                            double lr[] = {R(0) - dx, R(1) - 0.5 * dy};
                            double ur[] = {R(0)     , R(1) + 0.5 * dy};
                            *reinterpret_cast<double*>(&cparams[0] + 8 * sizeof(double)) = lr[0];
                            *reinterpret_cast<double*>(&cparams[0] + 9 * sizeof(double)) = ur[0];

                            eyVel.params = &cparams[0];

                            double abserr = 0.0;
                            gsl_integration_qag(&eyVel,
                                                lr[1],
                                                ur[1],
                                                1e-6,
                                                1e-3,
                                                LIMIT_INTEGRATION_SIZE,
                                                N_INTEG_SAMPLES,
                                                w,
                                                &intoneYInt,
                                                &abserr);
                        }
                        velAFD(i,j) += chargeSign * Vector_t(intoneXInt, intoneYInt);
                    }
                }
            }
        }

        Communicator::communicateFields(velAFD, localAFDomains, localEFDomains, true);

        for (int jo = 0; jo <= Ny; ++ jo) {
            Vector_t Pos(0.0, (jo - Ny / 2) * dy);
            double chargeY = -0.5e-9 * (std::erf((Pos(1) + 0.5 * dy) * convFactorInv) -
                                        std::erf((Pos(1) - 0.5 * dy) * convFactorInv));

            for (int io = 0; io <= Nx; ++ io) {
                Pos(0) = (io - Nx / 2) * dx;

                charge = (0.5 * chargeSign *
                          (std::erf((Pos(0) + 0.5 * dx) * convFactorInv) -
                           std::erf((Pos(0) - 0.5 * dx) * convFactorInv)) * chargeY * tcharge.getValue());


                for (int j = lDom[1].first(); j <= lDom[1].last(); ++ j) {
                    elem[1] = Index(j,j);

                    for (int i = lDom[0].first(); i <= lDom[0].last(); ++ i) {
                        elem[0] = Index(i,i);

                        velEFD.localElement(elem) += charge * velPrefactor * velAFD(i - io + Nx, j - jo + Ny);
                    }
                }
            }
        }
    }
    msg << "processing... done" << endl;

    gsl_integration_workspace_free(w);
    gsl_integration_workspace_free(v);

    IpplTimings::stopTimer(crunchTimer);

    driftEFD.fillGuardCells();
    velEFD.fillGuardCells();
    accEFD.fillGuardCells();

    std::cout.flags(ff);
    std::cout.precision(ss);

    Tenzor<double, DIM> M(0.0);
    M(0,0) = cos(phi);
    M(0,1) = sin(phi);
    M(1,0) = -sin(phi);
    M(1,1) = cos(phi);

    std::vector<double> colSPos(2 * Nx - 1, 0.0);
    std::vector<double> colField(6 * (2 * Nx - 1), 0.0);

    Vector_t dir(cos(phi), sin(phi));
    std::vector<Vector_t> scatterHelper;
    ParticleAttrib<Vector_t> scatterPosition;
    ParticleAttrib<Vector_t> scatterAccField;
    ParticleAttrib<Vector_t> scatterVelField;
    ParticleAttrib<Vector_t> scatterDriftField;
    int firstNumSample[] = {-1, 0};

    for (int i = 0; i < 2 * Nx - 1; ++ i) {
        Vector_t Pos = (i - Nx + 1.0) * dx * dir + bunchPosition;
        colSPos[i] = (i - Nx + 1.0) * dx;

        elem = velEFD.get_mesh().getCellContaining(Pos);
        if (lDom.contains(elem)) {
            if (firstNumSample[0] == -1) firstNumSample[0] = i;
            scatterHelper.push_back(Pos);
        }
    }
    firstNumSample[1] = scatterHelper.size();

    if (scatterHelper.size() > 0) {
        scatterPosition.create(scatterHelper.size());
        scatterAccField.create(scatterHelper.size());
        scatterVelField.create(scatterHelper.size());
        scatterDriftField.create(scatterHelper.size());
    }

    Vector_t centerPosition(bunchPosition);
    if (centeredVTK.getValue() && VTKoutput.getValue()) {
        origin -= bunchPosition;
        driftEFD.get_mesh().set_origin(origin);
        velEFD.get_mesh().set_origin(origin);
        accEFD.get_mesh().set_origin(origin);
        centerPosition = 0.0;
    }

    for (unsigned int i = 0; i < scatterHelper.size(); ++ i) {
        scatterPosition[i] = scatterHelper[i] - bunchPosition + centerPosition;
        scatterAccField[i] = 0.0;
        scatterVelField[i] = 0.0;
        scatterDriftField[i] = 0.0;
    }

    scatterAccField.gather(accEFD, scatterPosition, IntOp_t());
    scatterVelField.gather(velEFD, scatterPosition, IntOp_t());
    scatterDriftField.gather(driftEFD, scatterPosition, IntOp_t());

    for (unsigned int i = 0; i < scatterHelper.size(); ++ i) {
        scatterAccField[i] = dot(M, scatterAccField[i]);
        colField[6 * (i + firstNumSample[0])] = scatterAccField[i](0);
        colField[6 * (i + firstNumSample[0]) + 1] = scatterAccField[i](1);

        scatterVelField[i] = dot(M, scatterVelField[i]);
        colField[6 * (i + firstNumSample[0]) + 2] = scatterVelField[i](0);
        colField[6 * (i + firstNumSample[0]) + 3] = scatterVelField[i](1);

        scatterDriftField[i] = dot(M, scatterDriftField[i]);
        colField[6 * (i + firstNumSample[0]) + 4] = scatterDriftField[i](0);
        colField[6 * (i + firstNumSample[0]) + 5] = scatterDriftField[i](1);
    }

    int tag = Ippl::Comm->next_tag(IPPL_APP_TAG0, IPPL_APP_CYCLE);
    if (Ippl::myNode() == 0) {
        MPI_Status status;

        std::vector<int> firstNumSamples(2 * Ippl::getNodes(), 0);
        firstNumSamples[0] = firstNumSample[0];
        firstNumSamples[1] = firstNumSample[1];
        MPI_Gather(firstNumSample, 2, MPI_INT, &(firstNumSamples[0]), 2, MPI_INT, 0, MPI_COMM_WORLD);

        for (int i = 1; i < Ippl::getNodes(); ++ i) {
            if (firstNumSamples[2 * i + 1] > 0) {
                MPI_Recv(&(colField[6 * firstNumSamples[2 * i]]),
                         6 * firstNumSamples[2 * i + 1],
                         MPI_DOUBLE,
                         i,
                         tag,
                         MPI_COMM_WORLD,
                         &status);
            }
        }

        std::stringstream fname;
        fname << "FieldOnLine"
              << "_phi=" << std::fixed << std::setprecision(2) << phiBend.getValue()
              << "_L=" << std::fixed << std::setprecision(2) << L << ".dat";

        std::ofstream dataOut(fname.str());
        for (int i = 0; i < 2 * Nx - 1; ++ i) {
            dataOut << std::setw(18) << std::setprecision(8) << colSPos[i]
                    << std::setw(18) << std::setprecision(8) << colField[6 * i]
                    << std::setw(18) << std::setprecision(8) << colField[6 * i + 1]
                    << std::setw(18) << std::setprecision(8) << colField[6 * i + 2]
                    << std::setw(18) << std::setprecision(8) << colField[6 * i + 3]
                    << std::setw(18) << std::setprecision(8) << colField[6 * i + 4]
                    << std::setw(18) << std::setprecision(8) << colField[6 * i + 5]
                    << std::setw(18) << std::setprecision(8) << (-0.25e-9 * (std::erf((colSPos[i] + 0.5 * dx) * convFactorInv) -
                                                                             std::erf((colSPos[i] - 0.5 * dx) * convFactorInv)) * tcharge.getValue() * (std::erf(0.5 * dy * convFactorInv) - std::erf(-0.5 * dy * convFactorInv)))

                    << std::endl;
        }
        dataOut.close();

    } else {
        MPI_Gather(firstNumSample, 2, MPI_INT, NULL, 2, MPI_INT, 0, MPI_COMM_WORLD);
        if (firstNumSample[1] > 0) {
            MPI_Send(&(colField[6 * firstNumSample[0]]),
                     6 * firstNumSample[1],
                     MPI_DOUBLE,
                     0,
                     tag,
                     MPI_COMM_WORLD);
        }
    }



    if (VTKoutput.getValue()) {
        BinaryVtkFile result;
        result.addVectorField(driftEFD, "Edrift-Field");
        result.addVectorField(velEFD, "Evel-Field");
        result.addVectorField(accEFD, "Eacc-Field");
        driftEFD[lDom[0]][lDom[1]] += velEFD[lDom[0]][lDom[1]] + accEFD[lDom[0]][lDom[1]];
        result.addVectorField(driftEFD, "E-Field");
        std::stringstream fname;
        fname << "result_phi=" << std::fixed << std::setprecision(2) << phiBend.getValue()
              << "_L=" << L;

        result.writeFile(fname.str());
    }

    IpplTimings::stopTimer(mainTimer);
    IpplTimings::print();

    return 0;
}

double intExBendAcc(double x,
                    void* params)
{
    double* doubleParams = (double*) params;
    char* charParams = (char*) params;
    doubleParams[0] = x;
    double ly = doubleParams[8];
    double uy = doubleParams[9];
    gsl_integration_workspace* w = *reinterpret_cast<gsl_integration_workspace**>(charParams + 10 * sizeof(double));

    gsl_function exAcc;
    exAcc.function = &ExBendAcc;
    exAcc.params = params;

    double res = 0.0;
    double abserr = 0.0;
    gsl_integration_qag(&exAcc,
                        ly,
                        uy,
                        1e-4,
                        1e-3,
                        LIMIT_INTEGRATION_SIZE,
                        N_INTEG_SAMPLES,
                        w,
                        &res,
                        &abserr);

    return res;
}

double intEyBendAcc(double y,
                    void* params)
{
    double* doubleParams = (double*) params;
    char* charParams = (char*) params;
    doubleParams[0] = y;
    double lx = doubleParams[8];
    double ux = doubleParams[9];
    gsl_integration_workspace* w = *reinterpret_cast<gsl_integration_workspace**>(charParams + 10 * sizeof(double));

    gsl_function eyAcc;
    eyAcc.function = &EyBendAcc;
    eyAcc.params = params;

    double res = 0.0;
    double abserr = 0.0;
    gsl_integration_qag(&eyAcc,
                        lx,
                        ux,
                        1e-3,
                        1e-3,
                        LIMIT_INTEGRATION_SIZE,
                        N_INTEG_SAMPLES,
                        w,
                        &res,
                        &abserr);

    return res;
}

double intExBendVel(double x,
                    void* params)
{
    double * doubleParams = (double*) params;
    char* charParams = (char*) params;
    doubleParams[0] = x;
    double ly = doubleParams[8];
    double uy = doubleParams[9];
    gsl_integration_workspace* w = *reinterpret_cast<gsl_integration_workspace**>(charParams + 10 * sizeof(double));

    gsl_function exVel;
    exVel.function = &ExBendVel;
    exVel.params = params;

    double res = 0.0;
    double abserr = 0.0;
    gsl_integration_qag(&exVel,
                        ly,
                        uy,
                        1e-6,
                        1e-3,
                        LIMIT_INTEGRATION_SIZE,
                        N_INTEG_SAMPLES,
                        w,
                        &res,
                        &abserr);

    return res;
}

double intEyBendVel(double y,
                    void* params)
{
    double* doubleParams = (double*) params;
    char* charParams = (char*) params;
    doubleParams[0] = y;
    double lx = doubleParams[8];
    double ux = doubleParams[9];
    gsl_integration_workspace* w = *reinterpret_cast<gsl_integration_workspace**>(charParams + 10 * sizeof(double));

    gsl_function eyVel;
    eyVel.function = &EyBendVel;
    eyVel.params = params;

    double res = 0.0;
    double abserr = 0.0;
    gsl_integration_qag(&eyVel,
                        lx,
                        ux,
                        1e-6,
                        1e-3,
                        LIMIT_INTEGRATION_SIZE,
                        N_INTEG_SAMPLES,
                        w,
                        &res,
                        &abserr);

    return res;
}

double EyBendAcc(double x, void* params)
{
    double * doubleParams = (double*) params;
    const Vector_t R(x, doubleParams[0]);
    const Vector_t V(doubleParams[1], doubleParams[2]);
    const Vector_t A(doubleParams[3], doubleParams[4]);
    const double gammasqr = doubleParams[5];
    const double mincdt = doubleParams[6];
    const double maxcdt = doubleParams[7];


    const double gamma = sqrt(gammasqr);
    const double dotRV = dot(R, V);
    const double gammasqrdotRV = gammasqr * dotRV;
    const double dotRR = dot(R, R);
    const double asqr = gammasqrdotRV*dotRV + dotRR;
    const double cdt = std::max(mincdt, gammasqrdotRV + sqrt(gammasqr * asqr));
    if (cdt > maxcdt) return 0.0;

    const double aInv = 1.0 / sqrt(asqr);
    const double gammasqrInv = 1.0 / gammasqr;
    const double zsqr = cdt*cdt * gammasqrInv - (2 * dotRV * cdt + dotRR);
    const double zmaxsqr = maxcdt*maxcdt * gammasqrInv - (2 * dotRV * maxcdt + dotRR);
    const double z = sqrt(zsqr);
    const double zmax = sqrt(zmaxsqr);
    const double sqrtasqrzsqrInv = 1.0 / sqrt(asqr + zsqr);
    const double sqrtasqrzmaxsqrInv = 1.0 / sqrt(asqr + zmaxsqr);
    const double RcrossBdot_z = R(0) * A(1) - R(1) * A(0);

    if (gammasqr * zsqr < 1e-12) {
        const double B = -zmax * sqrtasqrzmaxsqrInv + asinh(zmax * aInv);
        const double C = atan(zmax * aInv) * aInv;
        const double D = zmax * sqrtasqrzmaxsqrInv * aInv*aInv;

        const double retVal = 2.0 * (-A(1) * B -
                                     gamma * V(0) * RcrossBdot_z * C -
                                     (gammasqrdotRV * V(0) + R(0)) * RcrossBdot_z * D) / gamma;

        return retVal;
    } else {
        const double B = -zmax * sqrtasqrzmaxsqrInv + asinh(zmax * aInv) + z * sqrtasqrzsqrInv - asinh(z * aInv);
        const double C = (atan(zmax * aInv) - atan(z * aInv)) * aInv;
        const double D = (zmax * sqrtasqrzmaxsqrInv - z * sqrtasqrzsqrInv) * aInv*aInv;

        const double retVal = 2.0 * (-A(1) * B -
                                     gamma * V(0) * RcrossBdot_z * C -
                                     (gammasqrdotRV * V(0) + R(0)) * RcrossBdot_z * D) / gamma;

        return retVal;
    }
}

double EyBendVel(double x, void* params)
{
    double * doubleParams = (double*) params;
    const Vector_t R(x, doubleParams[0]);
    const Vector_t V(doubleParams[1], doubleParams[2]);
    const double gammasqr = doubleParams[5];
    const double mincdt = doubleParams[6];
    const double maxcdt = doubleParams[7];

    const double dotRV = dot(R, V);
    const double gammasqrdotRV = gammasqr * dotRV;
    const double dotRR = dot(R, R);
    const double asqr = gammasqrdotRV*dotRV + dotRR;
    const double cdt = std::max(mincdt, gammasqrdotRV + sqrt(gammasqr * asqr));
    if (cdt >= maxcdt) return 0.0;

    const double gammasqrzsqr = cdt*cdt - gammasqr * (2 * dotRV * cdt + dotRR);
    const double gammasqrzmaxsqr = maxcdt*maxcdt - gammasqr * (2 * dotRV * maxcdt + dotRR);

    if (gammasqrzsqr < 1e-12) {
        return 2.0 * R(1) * sqrt(gammasqr * gammasqrzmaxsqr / (asqr*asqr * (gammasqr * asqr + gammasqrzmaxsqr)));
    } else {
        return 2.0 * R(1) * (sqrt(gammasqr * gammasqrzmaxsqr / (asqr*asqr * (gammasqr * asqr + gammasqrzmaxsqr))) -
                             sqrt(gammasqr * gammasqrzsqr    / (asqr*asqr * (gammasqr * asqr + gammasqrzsqr))));
    }
}

double ExBendAcc(double y, void* params)
{
    double * doubleParams = (double*) params;
    const Vector_t R(doubleParams[0], y);
    const Vector_t V(doubleParams[1], doubleParams[2]);
    const Vector_t A(doubleParams[3], doubleParams[4]);
    const double gammasqr = doubleParams[5];
    const double mincdt = doubleParams[6];
    const double maxcdt = doubleParams[7];

    const double gamma = sqrt(gammasqr);
    const double dotRV = dot(R, V);
    const double gammasqrdotRV = gammasqr * dotRV;
    const double dotRR = dot(R, R);
    const double asqr = gammasqrdotRV*dotRV + dotRR;
    const double cdt = std::max(mincdt, gammasqrdotRV + sqrt(gammasqr * asqr));
    if (cdt > maxcdt) return 0.0;

    const double aInv = 1.0 / sqrt(asqr);
    const double gammasqrInv = 1.0 / gammasqr;
    const double zsqr = cdt*cdt * gammasqrInv - (2 * dotRV * cdt + dotRR);
    const double zmaxsqr = maxcdt*maxcdt * gammasqrInv - (2 * dotRV * maxcdt + dotRR);
    const double z = sqrt(zsqr);
    const double zmax = sqrt(zmaxsqr);
    const double sqrtasqrzsqrInv = 1.0 / sqrt(asqr + zsqr);
    const double sqrtasqrzmaxsqrInv = 1.0 / sqrt(asqr + zmaxsqr);
    const double RcrossBdot_z = R(0) * A(1) - R(1) * A(0);

    if (gammasqr * zsqr < 1e-12) {
        const double B = -zmax * sqrtasqrzmaxsqrInv + asinh(zmax * aInv);
        const double C = atan(zmax * aInv) * aInv;
        const double D = zmax * sqrtasqrzmaxsqrInv * aInv*aInv;

        const double retVal = 2.0 * (-A(0) * B +
                                     gamma * V(1) * RcrossBdot_z * C +
                                     (gammasqrdotRV * V(1) + R(1)) * RcrossBdot_z * D) / gamma;

        return retVal;
    } else {
        const double B = -zmax * sqrtasqrzmaxsqrInv + asinh(zmax * aInv) + z * sqrtasqrzsqrInv - asinh(z * aInv);
        const double C = (atan(zmax * aInv) - atan(z * aInv)) * aInv;
        const double D = (zmax * sqrtasqrzmaxsqrInv - z * sqrtasqrzsqrInv) * aInv*aInv;

        const double retVal = 2.0 * (-A(0) * B +
                                     gamma * V(1) * RcrossBdot_z * C +
                                     (gammasqrdotRV * V(1) + R(1)) * RcrossBdot_z * D) / gamma;

        return retVal;
    }
}

double ExBendVel(double y, void* params)
{
    double * doubleParams = (double*) params;
    Vector_t R(doubleParams[0], y);
    Vector_t V(doubleParams[1], doubleParams[2]);
    double gammasqr = doubleParams[5];
    double mincdt = doubleParams[6];
    double maxcdt = doubleParams[7];

    const double dotRV = dot(R, V);
    const double gammasqrdotRV = gammasqr * dotRV;
    const double dotRR = dot(R, R);
    const double asqr = gammasqrdotRV*dotRV + dotRR;
    const double cdt = std::max(mincdt, gammasqrdotRV + sqrt(gammasqr * asqr));
    if (cdt >= maxcdt) return 0.0;

    const double gammasqrzsqr = cdt*cdt - gammasqr * (2 * dotRV * cdt + dotRR);
    const double gammasqrzmaxsqr = maxcdt*maxcdt - gammasqr * (2 * dotRV * maxcdt + dotRR);

    if (gammasqrzsqr < 1e-12) {
        return 2.0 * R(0) * sqrt(gammasqr * gammasqrzmaxsqr / (asqr*asqr * (gammasqr * asqr + gammasqrzmaxsqr)));
    } else {
        return 2.0 * R(0) * (sqrt(gammasqr * gammasqrzmaxsqr / (asqr*asqr * (gammasqr * asqr + gammasqrzmaxsqr))) -
                             sqrt(gammasqr * gammasqrzsqr    / (asqr*asqr * (gammasqr * asqr + gammasqrzsqr))));
    }
}


void errorHandler (const char * reason,
                   const char * file,
                   int line,
                   int gsl_errno)
{
    dbg << "ERROR " << file << ": " << line << "\t" << reason << "; " << gsl_errno << std::endl;
}
