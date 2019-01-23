//TODO: implement a dynamic ordinal numbering that will change depending what dofs are used .
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    /// <summary>
    /// Degrees of freedom corresponding to the displacement of a point/body along the 3 possible axes. Implements enum pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class DisplacementDof: IDof
    {
        /// <summary>
        /// This dof corresponds to the displacement along X axis.
        /// </summary>
        public static readonly DisplacementDof X = new DisplacementDof("Displacement along X axis");

        /// <summary>
        /// This dof corresponds to the displacement along Y axis.
        /// </summary>
        public static readonly DisplacementDof Y = new DisplacementDof("Displacement along Y axis");

        /// <summary>
        /// This dof corresponds to the displacement along Z axis.
        /// </summary>
        public static readonly DisplacementDof Z = new DisplacementDof("Displacement along Z axis");

        private readonly string name; 

        private DisplacementDof(string name) // No more displacement dofs can be created
        {
            this.name = name;
        }

        public override string ToString() => name;
    }
}
