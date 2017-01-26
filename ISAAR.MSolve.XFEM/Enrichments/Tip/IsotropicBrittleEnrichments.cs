using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Enrichments.Tip
{
    class IsotropicBrittleEnrichments
    {
        private class Psi1 : ITipEnrichment
        {
            public double ValueAt(double r, double a)
            {
                return Math.Sqrt(r) * Math.Sin(a / 2.0);
            }

            public double RadialDerivativeAt(double r, double a)
            {
                throw new NotImplementedException();
            }

            public double AngularDerivativeAt(double r, double a)
            {
                throw new NotImplementedException();
            }
        }

        private class Psi2 : ITipEnrichment
        {
            public double ValueAt(double r, double a)
            {
                return Math.Sqrt(r) * Math.Cos(a / 2.0);
            }

            public double RadialDerivativeAt(double r, double a)
            {
                throw new NotImplementedException();
            }

            public double AngularDerivativeAt(double r, double a)
            {
                throw new NotImplementedException();
            }
        }

        private class Psi3 : ITipEnrichment
        {
            public double ValueAt(double r, double a)
            {
                return Math.Sqrt(r) * Math.Sin(a / 2.0) * Math.Sin(a);
            }

            public double RadialDerivativeAt(double r, double a)
            {
                throw new NotImplementedException();
            }

            public double AngularDerivativeAt(double r, double a)
            {
                throw new NotImplementedException();
            }
        }

        private class Psi4 : ITipEnrichment
        {
            public double ValueAt(double r, double a)
            {
                return Math.Sqrt(r) * Math.Cos(a / 2.0) * Math.Sin(a);
            }

            public double RadialDerivativeAt(double r, double a)
            {
                throw new NotImplementedException();
            }

            public double AngularDerivativeAt(double r, double a)
            {
                throw new NotImplementedException();
            }
        }

        public static readonly ITipEnrichment First = new Psi1();
        public static readonly ITipEnrichment Second = new Psi2();
        public static readonly ITipEnrichment Third = new Psi3();
        public static readonly ITipEnrichment Fourth = new Psi4();

        private static ITipEnrichment[] AllFunctions
        {
            get { return new ITipEnrichment[] { First, Second, Third, Fourth }; }
        }

        public static double[] AllValuesAt(double r, double a)
        {
            double[] values = new double[4];
            ITipEnrichment[] functions = AllFunctions;
            for (int i = 0; i < 4; ++i)
            {
                values[i] = functions[i].ValueAt(r, a);
            }
            return values;
        }

        public static Tuple<double, double>[] AllDerivativesAt(double r, double a)
        {
            var derivatives = new Tuple<double, double>[4];
            ITipEnrichment[] functions = AllFunctions;
            for (int i = 0; i < 4; ++i)
            {
                derivatives[i] = new Tuple<double, double>(
                    functions[i].RadialDerivativeAt(r, a), functions[i].AngularDerivativeAt(r, a));
            }
            return derivatives;
        }
    }
}
