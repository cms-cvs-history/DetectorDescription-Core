#include <DetectorDescription/Core/src/EllipticalTube.h>
#include <DetectorDescription/Core/src/Sphere.h>
#include <DataFormats/GeometryVector/interface/Pi.h>
#include "CLHEP/Units/GlobalSystemOfUnits.h"
#include <G4EllipticalTube.hh>
#include <G4Sphere.hh>
#include <string>

void doEllipticalTube (const std::string& name, double dx, double dy, double dz) {

  std::cout << "elliptical tube major axis " << dx/cm << " cm; minor axis " << dy/cm << " cm; and " << dz/cm << " length" << std::endl;
  G4EllipticalTube g4t(name,dx, dy, dz);
  DDI::EllipticalTube ddt(dx, dy, dz);
  std::cout << "\tg4 volume = " << g4t.GetCubicVolume()/cm3 <<" cm3" << std::endl;
  std::cout << "\tdd volume = " << ddt.volume()/cm3 << " cm3"<<  std::endl;
  std::cout << "\tcalc volume = " << 2*dz*Geom::pi()*dy*dx / cm3 << " cm3 " <<std::endl;

}

void doSphere (const std::string& name, double rMin, double rMax, 
	       double pSPhi, double pDPhi, double pSTheta, double pDTheta ) {

  std::cout << "spherical shell rMin " << rMin/cm << " cm; rMax " << rMax/cm << " cm; ";
  std::cout << "startPhi " << pSPhi/deg << " deg; deltaPhi " << pDPhi/deg << " deg; "; 
  std::cout << "startTheta " << pSTheta/deg << " deg; deltaTheta " << pDTheta/deg << " deg; "; 
  std::cout << std::endl;
  G4Sphere g4(name,rMin, rMax, pSPhi, pDPhi, pSTheta, pDTheta);
  DDI::Sphere dd(rMin, rMax, pSPhi, pDPhi, pSTheta, pDTheta);
  std::cout << "\tg4 volume = " << g4.GetCubicVolume()/cm3 <<" cm3" << std::endl;
  std::cout << "\tdd volume = " << dd.volume()/cm3 << " cm3"<<  std::endl;

}

int main(int argc, char *argv[]) {

  std::cout << "\n\nElliptical Tube tests\n" << std::endl;
  double dx(2.*cm);
  double dy(2.*cm);
  double dz(2.*cm);
  std::string name("fred1");
  doEllipticalTube(name, dx, dy, dz);

  dy = 3.*cm;
  doEllipticalTube(name, dx, dy, dz);

  dx = 3.* cm;
  dy = 2.* cm;
  dz = 10.* cm;
  doEllipticalTube(name, dx, dy, dz);

  dx = 300.* cm;
  dy = 400.* cm;
  dz = 3000. * cm;
  doEllipticalTube(name, dx, dy, dz);

  std::cout << "\n\nSphere tests\n" << std::endl;
  std::cout << "\nThis next should be the same as a 2cm ball: " << std::endl;
  doSphere("fred1", 0.0*cm, 2.0*cm, 0.*deg, 360.*deg, 0., 180.*deg);
  std::cout << "Manual computation gives: " 
	    << 4./3. * Geom::pi() * 2.0*cm * 2.0*cm *2.0*cm / cm3
	    <<std::endl;
  std::cout << "If you mess up phi and theta you get: " << std::endl;
  doSphere("fred1", 0.0*cm, 2.0*cm, 0.*deg, 180.*deg, 0., 360.*deg);

  std::cout << "\n1 cm thick shell: " << std::endl;
  doSphere ("fred1", 2.0*cm, 3.0*cm, 0.*deg, 360.*deg, 0., 180.*deg);
  std::cout << "Manual computation gives: "
	    << 4./3. * Geom::pi() * 3.0*cm * 3.0*cm *3.0*cm / cm3 - 4./3. * Geom::pi() * 2.0*cm * 2.0*cm *2.0*cm / cm3
	    <<std::endl;

  std::cout << "\nHALF of the above 1 cm thick shell: " << std::endl;
  doSphere ("fred1", 2.0*cm, 3.0*cm, 0.*deg, 180.*deg, 0., 180.*deg);
  std::cout << "Manual computation gives: "
	    << (4./3. * Geom::pi() * 3.0*cm * 3.0*cm *3.0*cm / cm3 - 4./3. * Geom::pi() * 2.0*cm * 2.0*cm *2.0*cm / cm3) / 2.
	    <<std::endl;



  return EXIT_SUCCESS;
}

