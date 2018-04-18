using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Integration.Points;

namespace ISAAR.MSolve.XFEM.Materials
{
    static class MaterialUtilities
    {
        public static void CheckYoungModulus(double youngModulus)
        {
            if (youngModulus <= 0.0)
            {
                throw new ArgumentException("Young's modulus must be positive but was: " + youngModulus);
            }
        }

        public static void CheckPoissonRatio(double poissonRatio)
        {
            if (poissonRatio < -1.0 || poissonRatio > 0.5)
            {
                throw new ArgumentException("Poisson's ratio must be in the range [-1, 0.5] but was: " + poissonRatio);
            }
        }

        public static void CheckThickness(double thickness)
        {
            if (thickness <= 0.0) throw new ArgumentException("Thickness must be positive but was: " + thickness);
        }
    }
}
