using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using ISAAR.MSolve.Discretization.Interfaces;

namespace ISAAR.MSolve.FEM.Entities
{
    public class ElementMassAccelerationLoad_v2
    {
        public Element_v2 Element { get; set; }
        public DOFType DOF { get; set; }
        public double Amount { get; set; }
    }
}
