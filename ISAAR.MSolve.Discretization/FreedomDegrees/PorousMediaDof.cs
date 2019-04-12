namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    /// <summary>
    /// Degree of freedom corresponding to the pore pressure at a single point. Implements enum pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class PorousMediaDof : IDofType
    {
        /// <summary>
        /// This dof corresponds to the temperature.
        /// </summary>
        public static readonly PorousMediaDof Pressure = new PorousMediaDof("Pressure");

        private readonly string name;

        private PorousMediaDof(string name) => this.name = name;

        public override string ToString() => name;
    }
}
