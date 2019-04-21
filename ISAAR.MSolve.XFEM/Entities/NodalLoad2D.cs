using ISAAR.MSolve.Discretization.FreedomDegrees;
using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Entities
{
    class NodalLoad2D
    {
        public Node Node { get; }
        public StructuralDof DofType { get; }
        public double Value { get; }

        public NodalLoad2D(Node node, StructuralDof dofType, double value)
        {
            this.Node = node;
            this.DofType = dofType;
            this.Value = value;
        }
    }
}
