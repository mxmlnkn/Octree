/*
 rm integratemain.exe; g++ integratemain.cpp -o integratemain.exe -Wall -Wextra -Wchar-subscripts -Wcomment -Wconversion -Wformat -Wmissing-braces -Wmissing-noreturn -Wparentheses -Wreturn-type -Wshadow -Wsign-compare -Wstrict-aliasing -Wuninitialized -Wunknown-pragmas -Wunreachable-code -Wno-unused-parameter -Wno-unused-variable -Werror -std=c++0x -I ../libs/ -g 2>&1 | grep 'error'; ./integratemain.exe
*/

#include <iostream>
#include <string>
#include <algorithm>
#include "math/Polynomial.h"

int main( void ) {
    std::string input, variables;
    std::cout << "With which variables do you want to operate? (E.g. default: 'x,y,z' if no input given, max. 3 allowed ):\n    ";
    std::getline(std::cin, variables);
    if ( variables.length() == 0 )
        variables = "x,y,z";
    std::cout << "\nPlease input a polynomial command, like Integrate[r+1,r,x,1]\n";
    std::cout << "  syntax: Integrate[<polynomial>,<integration variable>,<lower bound>,<upper bound>] where lower and upper bound can again be polynomials\n    ";
    std::getline(std::cin, input);

    size_t vars = std::count( variables.begin(), variables.end(), ',' ) + 1;
    //std::cout << "Counted " << vars << " variables\n";
    std::string magicvarname = ",fgnjuztgfhju";
    while ( vars < 3 ) {
        variables += magicvarname;
        magicvarname[ magicvarname.length()-1 ]++;
        ++vars;
    }
    //std::cout << "Variables: " << variables << " Command: " << input << "\n";
    std::cout << "\n => " << Pol3Vars( variables, input ) << "\n";

    return 0;
}

