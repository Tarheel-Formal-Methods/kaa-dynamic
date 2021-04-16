#include <iostream>
#include <string>
#include <exception>
#include <map>
#include <vector>

#include <kodiak.hpp>

using namespace std;
using namespace kodiak;

// persistant variables
static bool bernstein = true;
static int precision = -6;

static vector<Real> reals;
static map<string, int> varToRealIndex;

extern "C"
{
    // call once at initialization time
    void init()
    {
        Kodiak::init();
        Kodiak::set_safe_input(false);
    }

    void free_stack()
    {
      reals.clear();
      varToRealIndex.clear();
    }

    // use bernstein expansion for solving? (false = interval artithmetic)
    void useBernstein(int val)
    {
        bernstein = val ? true : false;
    }

    // set answer accuracte (10^prec), more negative = more accurate
    void setPrecision(int prec)
    {
        precision = prec;
    }

    // create a new variable, with a bigger index than before
    int addVariable(const char* name)
    {
        if (varToRealIndex.find(name) != varToRealIndex.end())
        {
            return 0;
        }

        reals.push_back(var(name));
        varToRealIndex[name] = reals.size() - 1;

        if (varToRealIndex.size() != reals.size())
        {
            throw runtime_error("addVariable() called after other Kodiak expressions were created. "
                                "All addVariable() calls MUST happen first.");
        }

        return 1;
    }

    // lookup a variable expression index by name (should have been prevoisly inserted with addVariable)
    int lookupVariable(const char* name)
    {
        auto i = varToRealIndex.find(name);

        if (i == varToRealIndex.end())
        {
            string msg = "lookupVariable() called with unknown variable: (did you call addVariable() first?): ";
            msg += name;
            throw runtime_error(msg);
        }

        return i->second;
    }

    // make a new expression for a double value
    int makeDouble(double d)
    {
        reals.push_back(val(approx(d)));

        return reals.size() - 1;
    }

    void checkIndex(const char* name, int i)
    {
        if (i < 0 || i > (int)reals.size())
        {
            char msg[256];
            snprintf(msg, sizeof(msg), "%s() called with out-of-bounds expression index: %d, valid range is [0, %d]",
                     name, i, (int)reals.size());

            throw runtime_error(msg);
        }
    }

    // make a new expression for a multiplication
    int makeMult(int a, int b)
    {
        checkIndex(__func__, a);
        checkIndex(__func__, b);

        reals.push_back(reals[a] * reals[b]);

        return reals.size() - 1;
    }

    // print the given expression to stdout
    void printExpression(int a)
    {
        checkIndex(__func__, a);

        cout << reals[a] << endl;
    }

    // make a new expression for a multiplication
    int makeAdd(int a, int b)
    {
        checkIndex(__func__, a);
        checkIndex(__func__, b);

        reals.push_back(reals[a] + reals[b]);

        return reals.size() - 1;
    }

    // make a new expression for a square
    int makeSq(int a)
    {
        checkIndex(__func__, a);

        reals.push_back(Sq(reals[a]));

        return reals.size() - 1;
    }

    // make a new expression for a sqrt
    int makeSqrt(int a)
    {
        checkIndex(__func__, a);

        reals.push_back(Sqrt(reals[a]));

        return reals.size() - 1;
    }

    // make a new expression for a pow (a^intp), note: intp is NOT an expression index
    int makeIntPow(int a, int intp)
    {
        checkIndex(__func__, a);

        reals.push_back(reals[a]^intp);

        return reals.size() - 1;
    }

    int makeSin(int a)
    {
        checkIndex(__func__, a);

        reals.push_back(Sin(reals[a]));

        return reals.size() - 1;
    }

    int makeCos(int a)
    {
        checkIndex(__func__, a);

        reals.push_back(Cos(reals[a]));

        return reals.size() - 1;
    }

    int makeAtan(int a)
    {
        checkIndex(__func__, a);

        reals.push_back(Atan(reals[a]));

        return reals.size() - 1;
    }

    // [ctypes.c_int, # nonlinear expression
    // ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int,
    // ctypes.c_int, # bias
    // ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int, ctypes.c_int,
    // ndpointer(ctypes.c_double, flags="C_CONTIGUOUS"), ctypes.c_int]
    void minmax_diff(int nonlinearExp, double* linearApprox, int linearApproxSize,
                double bias, double* bounds, int boundsRows, int boundsCols, double* rv, int rvSize)
    {
        checkIndex(__func__, nonlinearExp);

        int numVars = (int)varToRealIndex.size();

        if (linearApproxSize != numVars)
        {
            char msg[256];
            snprintf(msg, sizeof(msg), "minmax_diff called with linearApprox with %d entries (system has %d variables)",
                     linearApproxSize, numVars);

            throw runtime_error(msg);
        }

        if (boundsRows != numVars || boundsCols != 2)
        {
            char msg[256];
            snprintf(msg, sizeof(msg), "minmax_diff() called with %dx%d bounds (expected %dx2)",
                     boundsRows, boundsCols, numVars);

            throw runtime_error(msg);
        }

        if (rvSize != 4)
        {
            char msg[256];
            snprintf(msg, sizeof(msg), "minmax_diff() called with rv array of size %d (expected 4)",
                     rvSize);

            throw runtime_error(msg);
        }

        /////////////////////////////
        Real linear = val(approx(bias));

        // variables are first n values in real vector
        for (int i = 0; i < numVars; ++i)
            linear = linear + val(approx(linearApprox[i])) * reals[i];

        Real diff = reals[nonlinearExp] - linear;

        MinMaxSystem sys;

        //cout << "Setting default enclosure method. \n";
        //cout << std::flush;

        sys.setDefaultEnclosureMethodTrueForBernsteinAndFalseForIntervalArithmetic(bernstein);
        sys.set_precision(precision);

        // assign bounds
        for (auto it = varToRealIndex.begin(); it != varToRealIndex.end(); ++it)
        {
            int row = it->second;

            // sys.var("x", 1.25, 1.55);
            const char* var = it->first.c_str();

            sys.var(var, approx(bounds[row * 2]), approx(bounds[row * 2 + 1]));
        }

        //cout << "Calling sys.minmax \n";
        //cout << std::flush;

        sys.minmax(diff);

        MinMax answer = sys.answer();
        rv[0] = answer.lb_of_min();
        rv[1] = answer.ub_of_max();

        rv[2] = answer.ub_of_min();
        rv[3] = answer.lb_of_max();
    }
}
