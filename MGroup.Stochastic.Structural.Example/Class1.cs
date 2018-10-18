using System;

namespace MGroup.Stochastic.Structural.Example
{
    public class SolveProblem
    {
        public void Solve()
        {
            const int iterations = 1000;
            const double youngModulus = 2.1e8;

            var evaluator = new StructuralStochasticEvaluator(youngModulus);
            var m = new MonteCarlo(iterations, evaluator, evaluator);
            m.Evaluate();
        }
    }
}
