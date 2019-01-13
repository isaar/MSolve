using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

//TODO: Many boundary (including this) should depend on INode and IElement, in order to work for all discretization methods
namespace ISAAR.MSolve.FEM.Entities
{
    public class Load_v2
    {
        public Node_v2 Node { get; set; }
        public DOFType DOF { get; set; }
        public double Amount { get; set; }
    }
}
