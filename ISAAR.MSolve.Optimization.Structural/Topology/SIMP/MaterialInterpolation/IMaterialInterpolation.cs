using System;
using System.Collections.Generic;
using System.Text;

//TODO: what about the young modulus defined in the element materials?
namespace ISAAR.MSolve.Optimization.Structural.Topology.SIMP.MaterialInterpolation
{
    public interface IMaterialInterpolation
    {
        /// <summary>
        /// The maximum allowable density of an element.
        /// </summary>
        double MaxDensity { get; }

        /// <summary>
        /// The minimum allowable density of an element.
        /// </summary>
        double MinDensity { get; }

        /// <summary>
        /// Calculates the interpolated material property of an element with density = <paramref name="elementDensity"/>.
        /// </summary>
        /// <param name="elementDensity">The element's density after any filter has been applied.</param>
        double CalcMaterialProperty(double elementDensity);

        /// <summary>
        /// Calculates the derivative of the interpolated material property of an element with respect to its density, when  
        /// density = <paramref name="elementDensity"/>.
        /// </summary>
        /// <param name="elementDensity">The element's density after any filter has been applied.</param>
        double CalcMaterialPropertyDerivative(double elementDensity);
    }
}
