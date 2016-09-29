using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISSAR.MSolve.IGAPreProcessor
{
    class ControlPoint
    {
        private int ID { get; }
        private double X { get; }
        private double Y { get; }
        private double Z { get; }
        private double Ksi { get; }
        private double Heta { get; }
        private double Zeta { get; }
        private double Weight { get; }

        public ControlPoint(int id, double x, double y, double z, double ksi, double heta, double zeta, double weight)
        {
            this.ID = id;
            this.X = x;
            this.Y = y;
            this.Z = z;
            this.Ksi = ksi;
            this.Heta = Heta;
            this.Zeta = zeta;
            this.Weight = weight;
        }
    }
}
