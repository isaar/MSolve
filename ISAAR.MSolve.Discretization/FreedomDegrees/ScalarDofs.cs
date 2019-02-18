//TODO: Perhaps each type should be a different class
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    /// <summary>
    /// Contains various freedom degrees that are not vectors.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public static class ScalarDofs
    {
        /// <summary>
        /// This dof corresponds to the pore pressure.
        /// </summary>
        public static readonly IDof PorePressure = new GenericScalarDof("Displacement along X axis");

        /// <summary>
        /// This dof corresponds to the temperature.
        /// </summary>
        public static readonly IDof Temperature = new GenericScalarDof("Displacement along Y axis");

        private class GenericScalarDof: IDof
        {
            private readonly string name;

            internal GenericScalarDof(string name)
            {
                this.name = name;
            }

            public override string ToString() => name;
        }
    }
}
