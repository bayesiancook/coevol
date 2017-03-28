#ifndef GENERICTIMELINE_H
#define GENERICTIMELINE_H

#include <iostream>
#include <map>
#include <vector>

class GenericTimeLine {
  public:
    GenericTimeLine() = default;
    virtual ~GenericTimeLine() = default;

    virtual double GetVal(double abstime, int TimeLineIndex) = 0;
    virtual int GetNumberOfTimeLines() = 0;

  protected:
    int NumberOfTimeLines;
};

/*****************************************************************/

class LinearlyInterpolatedTimeLine : public GenericTimeLine {
  public:
    LinearlyInterpolatedTimeLine() {}
    ~LinearlyInterpolatedTimeLine() override = default;

    double GetVal(double abstime, int TimeLineIndex) override {
        TimeLineIndex = TimeLineIndex - 1;
        if (TimeLineIndex > NumberOfTimeLines) {
            return (0);
        }
        std::map<double, double>::iterator Younger, Older;
        double ValueYounger, ValueOlder;
        Younger = TimesAndValues[TimeLineIndex].lower_bound(abstime);
        Older = TimesAndValues[TimeLineIndex].upper_bound(abstime);
        ValueYounger = (*Younger).second;
        if ((*Younger).first == abstime) {
            return ValueYounger;
        }
        if (Younger == TimesAndValues[TimeLineIndex].begin()) {
            return ValueYounger;
        } else if (Older == TimesAndValues[TimeLineIndex].end()) {
            Older--;
            return (*Older).second;
        } else if (Younger == Older) {
            Younger--;
            ValueYounger = (*(Younger)).second;
            ValueOlder = (*Older).second;
        } else if (abstime > (*TimesAndValues[TimeLineIndex].rbegin()).first) {
            return (*Older).second;
        } else {
            ValueOlder = (*Older).second;
        }
        return ValueYounger +
               (ValueOlder - ValueYounger) * (abstime - (*Younger).first) /
                   ((*Older).first - (*Younger).first);  // Linear interpolation
    }

    int GetNumberOfTimeLines() override { return NumberOfTimeLines; }

    // The format is going to be:
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
    void ToStream(std::ostream &os) const {
        std::map<double, double>::iterator it;
        std::map<double, double> temp;
        os << NumberOfTimeLines << std::endl;
        for (const auto &TimesAndValue : TimesAndValues) {
            temp = TimesAndValue;
            os << temp.size() << std::endl;
            for (it = temp.begin(); it != temp.end(); it++) {
                os << (*it).first << "\t" << (*it).second << std::endl;
            }
        }
        return;
    }

    void FromStream(std::istream &is) {
        double Age;
        double Value;
        int NumberOfValues;
        std::pair<double, double> ageValue;
        std::map<double, double> TimeValueMap;
        is >> NumberOfTimeLines;
        for (int i = 0; i < NumberOfTimeLines; i++) {
            is >> NumberOfValues;
            TimesAndValues.push_back(TimeValueMap);
            for (int j = 0; j < NumberOfValues; j++) {
                is >> Age >> Value;
                TimesAndValues[i][Age] = Value;
            }
        }
        return;
    }

  private:
    std::vector<std::map<double, double>> TimesAndValues;
};

#endif
