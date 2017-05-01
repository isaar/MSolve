using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Jintegral
{
    // TODO 1: An IHomogeneousMaterialField would be better than the restrictive HomogeneousElasticMaterial2D
    // TODO 2: Enforce that all elements of the (J-integral) domain have the identical material properties to the ones 
    //      passed to this class.
    class HomogeneousSIFCalculator : ISIFCalculator
    {
        private readonly double equivalentYoungModulus;

        /// <summary>
        /// The material properties (E, v, E*, v*) must be the same across all elements. The user assumes 
        /// responsibility for passing a <see cref="HomogeneousElasticMaterial2D"/> that has the same properties as 
        /// the materials of all other elements of the integration domain.
        /// </summary>
        /// <param name="globalMaterial">The material properties which must be identical for all elements and this 
        ///     class</param>
        public HomogeneousSIFCalculator(HomogeneousElasticMaterial2D globalMaterial)
        {
            this.equivalentYoungModulus = globalMaterial.HomogeneousEquivalentYoungModulus;
        }

        public double CalculateSIF(double interactionIntegral)
        {
            return 0.5 * equivalentYoungModulus * interactionIntegral;
        }
    }
}
