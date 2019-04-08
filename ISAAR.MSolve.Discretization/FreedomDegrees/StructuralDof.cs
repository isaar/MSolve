//TODO: implement a dynamic ordinal numbering that will change depending what dofs are used .
namespace ISAAR.MSolve.Discretization.FreedomDegrees
{
    /// <summary>
    /// Degrees of freedom corresponding to the movement of a point/body along or around the 3 possible axes. Implements enum 
    /// pattern.
    /// Authors: Serafeim Bakalakos
    /// </summary>
    public class StructuralDof: IDofType
    {
        /// <summary>
        /// This dof corresponds to translation along X axis.
        /// </summary>
        public static readonly StructuralDof TranslationX = new StructuralDof("TranslationX.");

        /// <summary>
        /// This dof corresponds to translation along Y axis.
        /// </summary>
        public static readonly StructuralDof TranslationY = new StructuralDof("TranslationY");

        /// <summary>
        /// This dof corresponds to translation along Z axis.
        /// </summary>
        public static readonly StructuralDof TranslationZ = new StructuralDof("TranslationZ");

        /// <summary>
        /// This dof corresponds to rotation around X axis.
        /// </summary>
        public static readonly StructuralDof RotationX = new StructuralDof("RotationX.");

        /// <summary>
        /// This dof corresponds to rotation around Y axis.
        /// </summary>
        public static readonly StructuralDof RotationY = new StructuralDof("RotationY");

        /// <summary>
        /// This dof corresponds to rotation around Z axis.
        /// </summary>
        public static readonly StructuralDof RotationZ = new StructuralDof("RotationZ");

        private readonly string name; 

        private StructuralDof(string name) => this.name = name;

        public override string ToString() => name;
    }
}
