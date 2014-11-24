package plugins.adufour.activecontours;

import icy.math.ArrayMath;

/**
 * Utility class defining a fixed-size window where a user may store values and check convergence
 * against various criteria
 * 
 * @author Alexandre Dufour
 */
public class SlidingWindow
{
    /**
     * The list of operations that can be applied on the window
     * 
     * @author Alexandre Dufour
     */
    public enum Operation
    {
        NONE, MIN, MAX, MEAN, SUM, VARIANCE,
        /**
         * Coefficient of variation (standard deviation over the mean)
         */
        VAR_COEFF
    };
    
    private double[] window;
    
    private int      count = 0;
    
    /**
     * Creates a new convergence window with given size, operation and convergence test sorting
     * method
     * 
     * @param size
     *            the window size
     */
    public SlidingWindow(int size)
    {
        setSize(size);
    }
    
    public int getSize()
    {
        return window.length;
    }
    
    public void setSize(int size)
    {
        window = new double[size];
        count = 0;
    }
    
    /**
     * Adds the given value to the queue
     * 
     * @param value
     */
    public final void push(double value)
    {
        // skip every other value to prevent oscillation effects
        if (count % 1 == 0) window[(count / 2) % window.length] = value;
        // window[count % window.length] = value;
        count++;
    }
    
    /**
     * Erase all values from the convergence window. Makes the window reusable without destruction
     */
    public void clear()
    {
        java.util.Arrays.fill(window, 0);
        count = 0;
    }
    
    public Double computeCriterion(Operation operation)
    {
        if (count < window.length * 2) return null;
        
        switch (operation)
        {
        case NONE:
            return null;
        case MIN:
            return ArrayMath.min(window);
        case MAX:
            return ArrayMath.max(window);
        case MEAN:
            return ArrayMath.mean(window);
        case SUM:
            return ArrayMath.sum(window);
        case VARIANCE:
            return ArrayMath.var(window, true);
        case VAR_COEFF:
            return ArrayMath.std(window, false) / ArrayMath.mean(window);
        default:
            throw new UnsupportedOperationException("operation " + operation.toString() + " not supported yet");
        }
    }
}
