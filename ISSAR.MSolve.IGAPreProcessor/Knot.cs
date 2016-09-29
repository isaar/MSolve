using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISSAR.MSolve.IGAPreProcessor
{
    class Knot
    {
        private int ID { get; }
        private int Ksi { get; }
        private int Heta { get; }
        private int Zeta { get; }

        public Knot(int id, int ksi, int heta, int zeta)
        {
            this.ID = id;
            this.Ksi = ksi;
            this.Heta = heta;
            this.Zeta = zeta;
        }
    }
}
