namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    /// <summary>
    /// Degree of freedom corresponding to the temperature at a single point. Implements enum pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class ThermalDof : IDofType
    {
        /// <summary>
        /// This dof corresponds to the temperature.
        /// </summary>
        public static readonly ThermalDof Temperature = new ThermalDof("Temperature");

        private readonly string name;

        private ThermalDof(string name) => this.name = name;

        public override string ToString() => name;
    }
}
