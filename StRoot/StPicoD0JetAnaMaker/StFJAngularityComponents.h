#ifndef ANGULARITY_COMPONENTS_H
#define ANGULARITY_COMPONENTS_H

#include "fastjet/PseudoJet.hh"
#include "fastjet/contrib/ShapeWithComponents.hh"

#include <vector>
#include <cmath>
#include <sstream>
#include <string>


class AngularityKappa1Components : public fastjet::contrib::ShapeWithComponents {
public:

  AngularityKappa1Components(double alpha = 1.0, double Radius = 0.4)
    : _alpha(alpha), _Radius(Radius)
  {
  }

  virtual std::string description() const
  {
    std::ostringstream oss;
    oss << "Generalized angularity lambda_alpha^1 with components, alpha = "
        << _alpha << ", R = " << _Radius;
    return oss.str();
  }

  virtual unsigned int n_components() const
  {
    return 2;
  }

  virtual std::vector<double> components(const fastjet::PseudoJet &jet) const
  {
    std::vector<double> comp(2, 0.0);

    double numerator = 0.0;

    std::vector<fastjet::PseudoJet> constituents = jet.constituents();

    for (size_t i = 0; i < constituents.size(); ++i) {

      // Same cut as in your current Angularity class:
      // keep charged + D0 + ghosts, skip unwanted inputs.
      if (constituents[i].user_index() < -1) continue;

      const double pt_i = constituents[i].pt();

      const double deta = jet.eta() - constituents[i].eta();
      double dphi = jet.phi() - constituents[i].phi();

      if (dphi >  M_PI) dphi -= 2.0 * M_PI;
      if (dphi < -M_PI) dphi += 2.0 * M_PI;

      const double dr = std::sqrt(deta*deta + dphi*dphi);

      numerator += pt_i * std::pow(dr / _Radius, _alpha);
    }

    // IMPORTANT:
    // To reproduce your current definition exactly, the denominator is jet.pt(),
    // because your original class uses:
    //
    //   pt_total = jet.pt()
    //
    // and not scalar sum_i pT_i.
    comp[0] = numerator;
    comp[1] = jet.pt();

    return comp;
  }

  virtual double result_from_components(const std::vector<double> &comp) const
  {
    const double numerator   = comp[0];
    const double denominator = comp[1];

    if (denominator <= 0.0) return -999.0;

    return numerator / denominator;
  }

private:
  double _alpha;
  double _Radius;
};

class MomentumDispersionComponents : public fastjet::contrib::ShapeWithComponents {
public:

  MomentumDispersionComponents()
  {
  }

  virtual std::string description() const
  {
    return "Momentum dispersion pTD with components";
  }

  virtual unsigned int n_components() const
  {
    return 2;
  }

  virtual std::vector<double> components(const fastjet::PseudoJet &jet) const
  {
    std::vector<double> comp(2, 0.0);

    double sumPt2 = 0.0;

    std::vector<fastjet::PseudoJet> constituents = jet.constituents();

    for (size_t i = 0; i < constituents.size(); ++i) {

      if (constituents[i].user_index() < -1) continue;

      const double pt_i = constituents[i].pt();

      sumPt2 += pt_i * pt_i;
    }

    comp[0] = sumPt2;
    comp[1] = jet.pt();

    return comp;
  }

  virtual double result_from_components(const std::vector<double> &comp) const
  {
    const double sumPt2 = comp[0];
    const double jetPt  = comp[1];

    if (sumPt2 < 0.0) return -999.0;
    if (jetPt <= 0.0) return -999.0;

    return std::sqrt(sumPt2) / jetPt;
  }

/*
virtual double result_from_components(const std::vector<double> &comp) const
{
  const double sumPt2 = comp[0];
  const double jetPt  = comp[1];

  if (!std::isfinite(sumPt2) || !std::isfinite(jetPt))
    return -999.0;

  // The normalized observable is undefined for a zero denominator.
  if (jetPt == 0.0)
    return -999.0;

  const double signedPtd =
      std::copysign(std::sqrt(std::abs(sumPt2)), sumPt2)
      / std::abs(jetPt);

  if (!std::isfinite(signedPtd))
    return -999.0;

  return signedPtd;
}
*/

};

#endif
