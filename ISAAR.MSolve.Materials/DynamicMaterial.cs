using System;
using System.Collections.Generic;
using System.Text;

//TODO: Density is also used in body forces, so it is not only for dynamic/modal problems.
//TODO: These cannot vary throughout the element. The solution is having a MaterialField that returns properties at specific 
//      natural points.
namespace ISAAR.MSolve.Materials
{
    /// <summary>
    /// Contains material properties used for dynamic or modal analysis. These are uniform for the whole element and immutable.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class DynamicMaterial
    {
        public DynamicMaterial(double density, double rayleighCoeffMass, double rayleighCoeffStiffness)
        {
            this.Density = density;
            this.RayleighCoeffMass = rayleighCoeffMass;
            this.RayleighCoeffStiffness = rayleighCoeffStiffness;
        }

        public double Density { get; }
        public double RayleighCoeffMass { get; }
        public double RayleighCoeffStiffness { get; }
    }
}
