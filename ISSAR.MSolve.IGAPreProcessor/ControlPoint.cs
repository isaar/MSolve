using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace ISSAR.MSolve.IGAPreProcessor
{
    public enum DOFType
    {
        Unknown = 0,
        X = 1,
        Y = 2,
        Z = 3
    }

    public class ControlPoint
    {
        public int ID { get; }
        public double X { get; }
        public double Y { get; }
        public double Z { get; }
        public double Ksi { get; }
        public double Heta { get; }
        public double Zeta { get; }
        public double Weight { get; }

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
