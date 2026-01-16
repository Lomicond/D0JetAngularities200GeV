#ifndef ANGULARITY_H
#define ANGULARITY_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/FunctionOfPseudoJet.hh"
#include <vector>
#include <cmath>
#include <sstream>
#include <string>

// Angularity class from fastjet::FunctionOfPseudoJet<double>
class Angularity : public fastjet::FunctionOfPseudoJet<double> {
public:
  //   a     - exponent for distance (delta_R)
  //   kappa - exponent for pT
  Angularity(double kappa = 1.0, double a = 1.0, double Radius = 0.4) :  _kappa(kappa), _a(a), _Radius(Radius) {}

  // Metoda, která spočítá generalized angularity pro daný jet
  virtual double result(const fastjet::PseudoJet& jet) const {
  
  double f_result = -999;
  
    double numerator = 0.0;
    double pt_total = jet.pt();
    std::vector<fastjet::PseudoJet> constituents = jet.constituents();
    
    for (size_t i = 0; i < constituents.size(); ++i) {
      double pt_i = constituents[i].pt();
      
      // Calculate the distance
      double deta = jet.eta() - constituents[i].eta();
      double dphi = jet.phi() - constituents[i].phi();
      
      // Angle range
      if (dphi > M_PI)  dphi -= 2 * M_PI;
      if (dphi < -M_PI) dphi += 2 * M_PI;
      
      double dr = std::sqrt(deta * deta + dphi * dphi);
      
      //Summing
      numerator += std::pow(1.*pt_i/pt_total, _kappa) * std::pow(1.0*dr/_Radius, _a);
    }
    
    
    f_result = numerator;
    
    return f_result;
    
  }

  virtual std::string description() const {
    std::ostringstream oss;
    oss << "Generalized Angularity with parameters kappa = " << _kappa << " and alpha = " << _a << " and jet R = " << _Radius;
    return oss.str();
  }
  
private:
  double _kappa; 
  double _a;     
  double _Radius;
};

#endif // ANGULARITY_H

