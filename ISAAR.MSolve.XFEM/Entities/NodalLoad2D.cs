using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using ISAAR.MSolve.XFEM.Entities.FreedomDegrees;

namespace ISAAR.MSolve.XFEM.Entities
{
    class NodalLoad2D
    {
        public Node2D Node { get; }
        public DisplacementDOF DofType { get; }
        public double Value { get; }

        public NodalLoad2D(Node2D node, DisplacementDOF dofType, double value)
        {
            this.Node = node;
            this.DofType = dofType;
            this.Value = value;
        }
    }
}
