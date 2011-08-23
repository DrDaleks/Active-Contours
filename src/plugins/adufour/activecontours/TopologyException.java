package plugins.adufour.activecontours;

/**
 * Class defining an exception which occurs when a contour is splitting during its evolution
 * 
 * @author Alexandre Dufour
 * 
 */
public class TopologyException extends Exception
{
	private static final long	serialVersionUID	= 1L;
	
	public final ActiveContour		source;
	
	public final ActiveContour[]		children;
	
	/**
	 * Creates a new Topology exception for the specified contour
	 * 
	 * @param contour
	 *            the contour undergoing a topology break
	 * @param children
	 *            an array containing zero or more contours that should replace the contour
	 *            raising the exception
	 */
	public TopologyException(ActiveContour contour, ActiveContour[] children)
	{
		super("Topology break detected in contour " + contour.hashCode());
		this.source = contour;
		this.children = children;
	}
	
}
