using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Quadratures;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    /// <summary>
    /// For linear problems with a single homogeneous material per element, 
    /// that does not change spatially (within the element) or temporally. Lightweight and fast.
    /// TODO: Needs a better name
    /// </summary>
    class SimpleIntegration2D: IIntegrationStrategy2D<ContinuumElement2D>
    {
        /// <summary>
        /// TODO: It must be immutable and its constitutive matrix must be fast to calculate: E.g. I could
        /// have a precached, immutable constitutive matrix. This would be inefficient if each element stored its 
        /// own constitutive matrix. Instead the constitutive matrix could be stored for the material instance and the
        /// <see cref="commonMaterial"/> fields in the different instances of <see cref="SimpleIntegration2D"/> (1 for each 
        /// element) could point to the same <see cref="IFiniteElementMaterial2D"/> object.
        /// 
        /// The above TODO doesn't take into account that materials store stresses, thus a different one is needed for
        /// each Gauss point.
        /// </summary>
        private readonly IFiniteElementMaterial2D commonMaterial;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="element"></param>
        /// <param name="commonMaterial">A homogeneous material, common for all gauss points. 
        ///     It will not change during the analysis</param>
        public SimpleIntegration2D(IFiniteElementMaterial2D commonMaterial)
        {
            this.commonMaterial = commonMaterial;
        }

        // TODO: Could be faster if the GaussPoint-Material pairs are stored instead, but that would require much more 
        // memory per element. Actually it would revert to the original design with GaussPoints to Materials 
        // dictionaries. Perhaps that dictionary cache would be faster for dynamic and quasi-static problems, but that 
        // could be done in another class! In static analysis the creation of that dictionary will be done exactly once,
        // so there is no need to cache it.
        public IReadOnlyDictionary<GaussPoint2D, IFiniteElementMaterial2D> GetIntegrationPointsAndMaterials(
            ContinuumElement2D element)
        {
            var pointsAndMaterials = new Dictionary<GaussPoint2D, IFiniteElementMaterial2D>();
            IReadOnlyList<GaussPoint2D> gaussPoints = element.StandardQuadrature.IntegrationPoints;
            foreach (GaussPoint2D gaussPoint in element.StandardQuadrature.IntegrationPoints)
            {
                pointsAndMaterials.Add(gaussPoint, commonMaterial.Clone());
            }
            return pointsAndMaterials;
        }

        public void Update(ContinuumElement2D element)
        {
            // Do nothing
        }
    }
}
