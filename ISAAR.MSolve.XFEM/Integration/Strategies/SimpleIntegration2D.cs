using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Elements;
using ISAAR.MSolve.XFEM.Integration.Points;
using ISAAR.MSolve.XFEM.Integration.Rules;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.Integration.Strategies
{
    /// <summary>
    /// For linear problems with a single homogeneous material per element, 
    /// that does not change spatially (within the element) or temporally. Lightweight and fast.
    /// TODO: Needs a better name
    /// </summary>
    class SimpleIntegration2D: IIntegrationStrategy2D
    {
        public class Factory: IIntegrationStrategyFactory2D
        {
            private readonly IFiniteElementMaterial2D material;

            public Factory(IFiniteElementMaterial2D material)
            {
                this.material = material;
            }

            public IIntegrationStrategy2D CreateStrategy(ContinuumElement2D elementType)
            {
                return new SimpleIntegration2D(elementType, material);
            }
        }

        /// <summary>
        /// TODO: It must be immutable and its constitutive matrix must be fast to calculate: E.g. I could
        /// have a precached, immutable constitutive matrix. This would be inefficient if each element stored its 
        /// own constitutive matrix. Instead the constitutive matrix could be stored for the material instance and the
        /// <see cref="material"/> fields in the different instances of <see cref="SimpleIntegration2D"/> (1 for each 
        /// element) could point to the same <see cref="IFiniteElementMaterial2D"/> object.
        /// </summary>
        private readonly IFiniteElementMaterial2D material;

        /// <summary>
        /// TODO: Create a standard integration rule interface that guarantees gauss points that are
        /// i) immutable, ii) precached for fast generation, iii) stored globally for all elements
        /// </summary>
        private readonly IIntegrationRule2D quadrature;

        /// <summary>
        /// 
        /// </summary>
        /// <param name="element"></param>
        /// <param name="material">A homogeneous material, common for all gauss points. 
        ///     It will not change during the analysis</param>
        private SimpleIntegration2D(ContinuumElement2D element, IFiniteElementMaterial2D material)
        {
            this.quadrature = element.StandardQuadrature;
            this.material = material;
        }

        // TODO: Could be faster if the GaussPoint-Material pairs are stored instead, but that would require much more 
        // memory per element. Actually it would revert to the original design with GaussPoints to Materials 
        // dictionaries. Perhaps that dictionary cache would be faster for dynamic and quasi-static problems, but that 
        // could be done in another class! In static classes the creation of that dictionary will be done exactly once,
        // so there is no need to cache it.
        public IEnumerable<Tuple<GaussPoint2D, IFiniteElementMaterial2D>> GetIntegrationPointsAndMaterials()
        {
            IReadOnlyList<GaussPoint2D> gaussPoints = quadrature.GenerateIntegrationPoints();
            var result = new Tuple<GaussPoint2D, IFiniteElementMaterial2D>[gaussPoints.Count];
            for (int i = 0; i < gaussPoints.Count; ++i)
            {
                result[i] = new Tuple<GaussPoint2D, IFiniteElementMaterial2D>(gaussPoints[i], material);
            }
            return result;
        }

        public void Update()
        {
            // Do nothing
        }
    }
}
