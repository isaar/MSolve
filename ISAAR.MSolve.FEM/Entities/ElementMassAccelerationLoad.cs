using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.FEM.Entities
{
    public class ElementMassAccelerationLoad
    {
        public Element Element { get; set; }
        public DOFType DOF { get; set; }
        public double Amount { get; set; }
    }
}
