using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISSAR.MSolve.IGAPreProcessor
{
    public class Knot
    {
        public int ID { get; }
        public double Ksi { get; }
        public double Heta { get; }
        public double Zeta { get; }

        public Knot(int id, double ksi, double heta, double zeta)
        {
            this.ID = id;
            this.Ksi = ksi;
            this.Heta = heta;
            this.Zeta = zeta;
        }
    }
}
