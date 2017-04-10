using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISAAR.MSolve.XFEM.CrackPropagation.Direction
{
    // TODO: Perhaps it should be (-pi/2, pi/2)
    interface ICrackGrowthDirectionLaw2D
    {
        /// <summary>
        /// Computes the angle coordinate of the next crack tip in the local polar system of the current crack tip.
        /// The angle belongs to the interval (-pi, pi]. However very sharp direction changes 
        /// are suspicious.
        /// </summary>
        /// <param name="sif1">The stress intensity factor of crack mode I (opening mode)</param>
        /// <param name="sif2">The stress intensity factor of crack mode II (sliding mode)</param>
        /// <returns>The angle belonging in the interval (-pi, pi] </returns>
        double ComputeGrowthAngle(double sifMode1, double sifMode2);
    }
}
