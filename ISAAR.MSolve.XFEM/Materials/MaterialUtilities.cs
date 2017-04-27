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
        public static IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> AssignMaterialToIntegrationPoints(
            IReadOnlyList<GaussPoint2D> integrationPoints, IFiniteElementMaterial2D commonMaterial)
        {
            var gpToMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            foreach (var point in integrationPoints)
            {
                gpToMaterials[point] = commonMaterial.Clone();
            }
            return gpToMaterials;
        }

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
