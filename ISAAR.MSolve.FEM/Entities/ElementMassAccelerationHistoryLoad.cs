using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;

namespace ISAAR.MSolve.FEM.Entities
{
    public class ElementMassAccelerationHistoryLoad
    {
        public Element Element { get; set; }
        public MassAccelerationHistoryLoad HistoryLoad { get; set; }
    }
}
