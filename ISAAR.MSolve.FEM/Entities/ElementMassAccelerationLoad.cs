using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.FreedomDegrees;

namespace ISAAR.MSolve.FEM.Entities
{
    public class ElementMassAccelerationLoad
    {
        public Element Element { get; set; }
        public IDofType DOF { get; set; }
        public double Amount { get; set; }
    }
}
