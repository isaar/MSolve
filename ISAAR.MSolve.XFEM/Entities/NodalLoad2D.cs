using ISAAR.MSolve.FEM.Entities;
using ISAAR.MSolve.XFEM.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Entities
{
    class NodalLoad2D
    {
        public Node Node { get; }
        public DisplacementDof DofType { get; }
        public double Value { get; }

        public NodalLoad2D(Node node, DisplacementDof dofType, double value)
        {
            this.Node = node;
            this.DofType = dofType;
            this.Value = value;
        }
    }
}
