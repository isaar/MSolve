namespace ISAAR.MSolve.Solvers.Commons
{
    public interface ISystemMatrixObserver
    {
        /// <summary>
        /// It will be called before setting the system matrix.
        /// </summary>
        void HandleMatrixWillBeSet();
    }
}
