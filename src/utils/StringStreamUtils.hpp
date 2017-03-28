#ifndef STRINGSTREAMUTILS_H
#define STRINGSTREAMUTILS_H

#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>

// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------
//     String Stream Utilities
// ---------------------------------------------------------------------------
// ---------------------------------------------------------------------------

const char digit[10] = {'0', '1', '2', '3', '4', '5', '6', '7', '8', '9'};

inline int GoPastNext(std::istream &is, const char inChar) {
    unsigned char c;
    if (!is.eof()) {
        do {
            is >> c;
        } while (c != inChar && !is.eof());
    }
    return static_cast<int>(!is.eof());
}

inline std::string ReadLine(std::istream &is) {
    std::string str = "";
    char c;
    do {
        is.get(c);
        if (c != '\n') {
            str += c;
        }
    } while (c != '\n' && !is.eof());
    return str;
}

inline void GoPastNextWord(std::istream &is, const std::string inWord) {
    unsigned int k = 0;
    char c;
    while ((!is.eof()) && (k < inWord.length())) {
        is.get(c);
        if ((c >= 65) && (c <= 90)) {
            c += 32;
        }
        char ca = inWord[k];
        if ((ca >= 65) && (ca <= 90)) {
            ca += 32;
        }
        if (c == ca) {
            k++;
        } else {
            k = 0;
        }
    }
}

inline int EquivalentStrings(std::string a, std::string b) {
    if (a.length() != b.length()) {
        return 0;
    }
    unsigned int k = 0;
    int cont = 1;
    while ((k < a.length()) && ((cont) != 0)) {
        char ca = a[k];
        char cb = b[k];
        if ((ca >= 65) && (ca <= 90)) {
            ca += 32;
        }
        if ((cb >= 65) && (cb <= 90)) {
            cb += 32;
        }
        if (ca != cb) {
            cont = 0;
        }
        k++;
    }
    return cont;
}

inline void GoPastNextLine(std::istream &is, const std::string inLine) {
    std::string theLine;
    do {
        theLine = ReadLine(is);
        std::cerr << theLine << "\n";
    } while (EquivalentStrings(theLine, inLine) == 0);
}

inline std::string StringReplace(char c, std::string by, std::string s) {
    std::string tmp;
    for (char i : s) {
        if (i == c) {
            tmp += by;
        } else {
            tmp += i;
        }
    }
    return tmp;
}

inline int EmptyLine(std::string s) {
    int unsigned n = 0;
    while ((n < s.length()) && ((s[n] == ' ') || (s[n] == '\t') || (s[n] == '\n'))) {
        n++;
    }
    return static_cast<int>(n == s.length());
}

inline std::string Filter(std::string input, char c) {
    std::string temp = "";
    for (char i : input) {
        if (i != c) {
            temp += i;
        }
    }
    return temp;
}

inline int IsInt(std::string s) {
    int returnValue = 1;
    unsigned int i = 0;
    if ((s[0] == '+') || (s[0] == '-')) {
        i++;
    }
    if (i == s.length()) {
        returnValue = 0;
    }

    while ((returnValue != 0) && (i < s.length())) {
        int j = 0;
        while ((j < 10) && (digit[j] != s[i])) {
            j++;
        }
        if (j == 10) {
            returnValue = 0;
        }
        i++;
    }
    return returnValue;
}

inline int IsFloat(std::string s) {
    int returnValue = 1;
    unsigned int i = 0;

    while ((returnValue != 0) && (i < s.length())) {
        int j = 0;
        while ((j < 10) && (digit[j] != s[i])) {
            j++;
        }
        if (j == 10) {
            if (!((s[i] == '.') || (s[i] == '-') || (s[i] == '+') || (s[i] == 'e'))) {
                returnValue = 0;
            }
        }
        i++;
    }
    return returnValue;
}

inline int IsDigit(char c) {
    int returnValue = 0;
    int i = 0;
    while ((returnValue == 0) && i < 10) {
        returnValue = static_cast<int>(c == digit[i]);
        i++;
    }
    return returnValue;
}

#endif
