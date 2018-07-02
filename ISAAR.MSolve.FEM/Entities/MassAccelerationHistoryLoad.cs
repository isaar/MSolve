using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.IO;
using System.Globalization;
using ISAAR.MSolve.Discretization.Interfaces;
using ISAAR.MSolve.FEM.Interfaces;

namespace ISAAR.MSolve.FEM.Entities
{
    public class MassAccelerationHistoryLoad : IMassAccelerationHistoryLoad
    {
        private readonly List<double> accelerationLoads = new List<double>();
        private readonly double magnifier = 1d;
        public virtual DOFType DOF { get; set; }

        public virtual double this[int currentTimeStep]
        {
            get { return currentTimeStep < accelerationLoads.Count ? accelerationLoads[currentTimeStep] * magnifier : 0; }
        }

        protected MassAccelerationHistoryLoad()
        {
        }

        public MassAccelerationHistoryLoad(string fileName, double magnifier)
        {
            this.magnifier = magnifier;
            using (StreamReader sr = new StreamReader(fileName))
                while (sr.Peek() >= 0)
                    accelerationLoads.Add(Double.Parse(sr.ReadLine(), new CultureInfo("en-US", false).NumberFormat));
        }

        public MassAccelerationHistoryLoad(string fileName) : this(fileName, 1d)
        {
        }
    }

    public class EmptyMassAccelerationHistoryLoad : MassAccelerationHistoryLoad
    {
        public EmptyMassAccelerationHistoryLoad() : base() { }
        public override DOFType DOF 
        {
            get { return DOFType.X; }
            set { } 
        }

        public override double this[int currentTimeStep] { get { return 0.0; } }
    }
}
