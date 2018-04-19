using System;
using System.Collections.Generic;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public class Constraint
    {
        public DOFType DOF { get; set; }//TODOMaria: make constraint immutable (add a constructor taking Node and DOF as arguments) (the same goes for Load, Node)

        public double Amount { get; set; }
    }
}

