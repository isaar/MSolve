using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public class ElementMassAccelerationHistoryLoad_v2
    {
        public Element_v2 Element { get; set; }
        public MassAccelerationHistoryLoad HistoryLoad { get; set; }
    }
}
