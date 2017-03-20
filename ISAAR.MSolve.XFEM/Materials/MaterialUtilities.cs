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
    }
}
