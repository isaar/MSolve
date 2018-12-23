namespace ISAAR.MSolve.Solvers.Commons
{
    public interface ISystemMatrixObserver
    {
        /// <summary>
        /// It will be called immediately before setting the system matrix.
        /// </summary>
        void OnMatrixSetting();
    }
}
