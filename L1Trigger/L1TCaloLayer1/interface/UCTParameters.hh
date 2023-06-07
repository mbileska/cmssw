#ifndef L1Trigger_L1TCaloLayer1_UCTParameters_hh
#define L1Trigger_L1TCaloLayer1_UCTParameters_hh

class UCTParameters {
public:

  UCTParameters(double activityFractionIn = 0.125,
		double ecalActivityFractionIn = 0.25,
		double miscActivityFractionIn = 0.25) :
    activityFractionValue(activityFractionIn),
    ecalActivityFractionValue(ecalActivityFractionIn),
    miscActivityFractionValue(miscActivityFractionIn)
  {;}

  virtual ~UCTParameters() {;}

  // More access functions

  const double activityFraction() {return activityFractionValue;}
  const double ecalActivityFraction() {return ecalActivityFractionValue;}
  const double miscActivityFraction() {return miscActivityFractionValue;}

  // Print parameters

  friend std::ostream& operator<<(std::ostream& os, const UCTParameters& p) {
    os << "UCTParmeters::activityFraction     = " << p.activityFractionValue     << std::endl;
    os << "UCTParmeters::ecalActivityFraction = " << p.ecalActivityFractionValue << std::endl;
    os << "UCTParmeters::miscActivityFraction = " << p.miscActivityFractionValue << std::endl;
    return os;
  }

private:

  // No copy constructor is needed

  UCTParameters(const UCTParameters&);

  // No equality operator is needed

  const UCTParameters& operator=(const UCTParameters&);

  // The parameters

  double activityFractionValue;
  double ecalActivityFractionValue;
  double miscActivityFractionValue;

};

#endif
