using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Materials;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Jintegral
{
    class HomogeneousSIFCalculator: ISIFCalculator
    {
        private readonly double equivalentYoungModulus;

        // TODO: when material fields are implemented, use a HomogeneousMaterialField class instead of a material point
        public HomogeneousSIFCalculator(IFiniteElementMaterial2D homogeneousMaterial)
        {
            this.equivalentYoungModulus = homogeneousMaterial.EquivalentYoungModulus;
        }

        public double CalculateSIF(double interactionIntegral)
        {
            return 0.5 * equivalentYoungModulus * interactionIntegral;
        }
    }
}
