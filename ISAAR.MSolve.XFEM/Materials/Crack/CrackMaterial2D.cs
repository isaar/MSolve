using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.Materials.Crack
{
    /// <summary>
    /// Perhaps It should also contain the SIFs. In the general case, a new material point is required at each time 
    /// step, along with the new SIFs.
    /// </summary>
    class CrackMaterial2D
    {
        public CrackMaterial2D(double poissonRatio, double shearModulus, double criticalFractureTougnhess)
        {
            // Arguments should be checked by the factory, to avoid rechecking at each step.

            this.PoissonRatio = poissonRatio;
            this.ShearModulus = shearModulus;
            this.CriticalFractureTougnhess = criticalFractureTougnhess;
        }

        double PoissonRatio { get; }
        double ShearModulus { get; }
        double CriticalFractureTougnhess { get; }
    }
}
