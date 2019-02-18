namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    /// <summary>
    /// Degrees of freedom corresponding to the rotation of a point/body around the 3 possible axes. Implements enum pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class RotationDof: IDof
    {
        /// <summary>
        /// This dof corresponds to the rotation around X axis.
        /// </summary>
        public static readonly RotationDof X = new RotationDof("Rotation around X axis");

        /// <summary>
        /// This dof corresponds to the rotation around Y axis.
        /// </summary>
        public static readonly RotationDof Y = new RotationDof("Rotation around Y axis");

        /// <summary>
        /// This dof corresponds to the rotation around Z axis.
        /// </summary>
        public static readonly RotationDof Z = new RotationDof("Rotation around Z axis");

        private readonly string name;

        private RotationDof(string name) // No more rotation dofs can be created
        {
            this.name = name;
        }

        public override string ToString() => name;
    }
}
