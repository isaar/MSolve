using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Termination
{
    class HomogeneousToughnessCriterion: IMaterialCriterion
    {
        private readonly double criticalFractureTougnhess;

        // TODO: when material fields are implemented, use a HomogeneousCrackMaterialField class to access the toughness
        public HomogeneousToughnessCriterion(double criticalFractureTougnhess)
        {
            // TODO: criticalFracture toughness should be checked by the corresponding material field
            this.criticalFractureTougnhess = criticalFractureTougnhess;
        }

        public bool Terminate(double sifMode1, double sifMode2)
        {
            double equivalentSIF = Math.Sqrt(sifMode1 * sifMode1 + sifMode2 * sifMode2);
            if (equivalentSIF < criticalFractureTougnhess) return false;
            else return true;
        }
    }
}
