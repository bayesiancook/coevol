#ifndef GENERICTIMELINE_H
#define GENERICTIMELINE_H

#include "StringStreamUtils.h"
#include <iostream>
#include <map>
#include <vector>
using namespace std;

class GenericTimeLine {
 
public:
  GenericTimeLine() {}
  virtual ~GenericTimeLine() {}
 
  virtual double GetVal(double abstime, int TimeLineIndex) = 0;
  virtual int GetNumberOfTimeLines() = 0;
 
protected:
  int NumberOfTimeLines;
 
 
};

/*****************************************************************/


class LinearlyInterpolatedTimeLine : public GenericTimeLine{
 
public:
  LinearlyInterpolatedTimeLine() {}
  virtual ~LinearlyInterpolatedTimeLine() {}
 
  double GetVal(double abstime, int TimeLineIndex) {
    TimeLineIndex = TimeLineIndex - 1;
    if (TimeLineIndex > NumberOfTimeLines) {
      return (0);
    }
    map < double,double>::iterator Younger, Older;
    double ValueYounger, ValueOlder;
    Younger = TimesAndValues[TimeLineIndex].lower_bound(abstime);
    Older = TimesAndValues[TimeLineIndex].upper_bound(abstime);
    ValueYounger = (*Younger).second;
    if ((*Younger).first == abstime) {
      return ValueYounger;
    }
    else if (Younger == TimesAndValues[TimeLineIndex].begin()) {
      return ValueYounger;
    }
    else if (Older == TimesAndValues[TimeLineIndex].end()) {
      Older--;
      return (*Older).second;
    }
    else if (Younger == Older) {
      Younger--;
      ValueYounger = (*(Younger)).second;
      ValueOlder = (*Older).second;
    }
    else if (abstime > (*TimesAndValues[TimeLineIndex].rbegin()).first) {
      return (*Older).second;
    }   
    else {
      ValueOlder = (*Older).second;
    }
    return  ValueYounger + (ValueOlder - ValueYounger) * ( abstime - (*Younger).first) / ((*Older).first - (*Younger).first) ; //Linear interpolation
  }
 
  int GetNumberOfTimeLines() {
    return NumberOfTimeLines;
  }

 
  //The format is going to be:
  /*
   numberOfTimelines
   numberOfValuesForTimeline1
   date1 value1
   date2 value2
   numberOfValuesForTimeline2
   date1 value1
   date2 value2
   ...
   */ 
  void ToStream(ostream& os) const {
    map < double, double>::iterator it;
    map < double, double> temp;
    os << NumberOfTimeLines << endl;
    for (unsigned int i = 0 ; i < TimesAndValues.size() ; i++) {
      temp = TimesAndValues[i];
      os <<  temp.size() <<endl;
      for ( it = temp.begin() ; it != temp.end(); it++ ) {
        os << (*it).first << "\t" << (*it).second <<  endl;
      }
    }
    return;
  }
 
	void FromStream(istream& is)	{
    double Age;
    double Value;
    int NumberOfValues;
    pair < double, double > ageValue;
    map<double, double> TimeValueMap;
        is >> NumberOfTimeLines;
        for (int i=0; i < NumberOfTimeLines; i++)	{     
          is >> NumberOfValues;
          TimesAndValues.push_back(TimeValueMap);
          for (int j = 0 ; j <  NumberOfValues ; j++ ) {
            is >> Age >> Value ;
            TimesAndValues[i][Age] = Value;
          }
        }
    return;
	}
 
  private:
  vector < map < double, double > > TimesAndValues;
 
 
 
};

#endif
