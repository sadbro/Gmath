class MetricError(Exception): pass

class Metrics:

    """
    ---------------------------------------------------------------------------
    This class is a implementation of a standard-ish Maths library.

                                                                    -- @sadbro
    ---------------------------------------------------------------------------
    Methods:

    Basic trigonometric functions like sin() and cos()
    Some constant generators like pi() and e()
    Some utility functions like fact() and sqrt()

    P.S. I know this is mostly incomplete and more functions will be added in the future :)

    ---------------------------------------------------------------------------

    """

    @staticmethod
    def pi(r= 6):

        iter=100000
        result= 0
        for i in range(iter):
            parity= (-1)**i
            deno=(2*i) +1
            result+= 4*(parity/deno)

        return round(result, r)

    @staticmethod
    def e(num=1, iter=10):

        result= 0
        for i in range(iter):
            parity= (2*i) +2
            deno= Metrics.fact((2*i) +1)
            result+= (parity/deno)

        return result**num

    @staticmethod
    def log(num, iter=10000, precision=3):

        if num > 0:
            N= (num-1)/(num+1)
            l= N
            for i in range(1, iter):
                K= (2*i)+1
                l+= (N**K)/K

            return round(2*l, precision)

        else:
            raise MetricError("Invalid parameter for log: non-positive entry `{}`".format(num))

    @staticmethod
    def prod(*num, s=1):

        for n in num:
            s*= n

        return s

    @staticmethod
    def sqrt(num, iter=1000, precision=3):

        x= 1
        for i in range(iter):
            x= 0.5*(x+(num/x))

        return round(x, precision)

    @staticmethod
    def fact(x):

        if x in [0, 1]:
            return 1

        else:
            result= x
            for i in range(1, x):
                result*= i

            return result

    @staticmethod
    def deg2rad(theta):

        return (Metrics.pi()/180)*theta

    @staticmethod
    def rad2deg(radian):

        return (180/Metrics.pi())*radian

    @staticmethod
    def sin(x, iter=80, points=4, unit="rad"):

        if unit == "rad":
            x_converted= x
        elif unit == "deg":
            x_converted= Metrics.deg2rad(x)
        else:
            raise TypeError("Unidentified angle unit: {}".format(unit))

        result= 0
        for i in range(iter):
            coef=(2*i) +1
            parity=(-1)**i
            result+= parity*((x_converted**coef)/Metrics.fact(coef))

        return round(result, points)

    @staticmethod
    def cos(x, iter=80, points=4, unit="rad"):

        if unit == "rad":
            x_converted= x
        elif unit == "deg":
            x_converted= Metrics.deg2rad(x)
        else:
            raise TypeError("Unidentified angle unit: {}".format(unit))

        result= 0
        for i in range(iter):
            coef=(2*i)
            parity=(-1)**i
            result+= parity*((x_converted**coef)/Metrics.fact(coef))

        return round(result, points)

    @staticmethod
    def countLR(array: list, num):

        fromL= 0
        for i in array:
            if i == num:
                fromL+= 1
                continue
            else:
                break

        fromR= 0
        for i in array[::-1]:
            if i == num:
                fromR+= 1
                continue
            else:
                break

        return fromL, fromR


    @staticmethod
    def anti_sparse(array: list):

        """
        returns a reduced array of sparse array {array}

        -------------------------------------------------------------------------
        Parameters: sparse array                 -> array

        -------------------------------------------------------------------------
        Examples:

        array= [0, 0, 1, 2 , 0, 0, 3, 1, 0, 7, 0, 0, 0, 0, 0, 0, 2, 0]
        anti_sparse(array) returns [(2, 1), (3, 2), (6, 3), (7, 1), (9, 7), (16, 2), (2, 1)], 0.67 <- sparsity

        This technique appends a tuple at the end containing the information about the no of 0s in
        the beginning and the end.

        -------------------------------------------------------------------------

        """

        sparsity= round(array.count(0)/len(array), 2)

        if type(array) != list:
            raise TypeError("Enter a valid type")
        else:
            res= []
            for i, e in enumerate(array):
                if e != 0:
                    res.append((i, e))

            res.append(Metrics.countLR(array, 0))

            return res, sparsity

    @staticmethod
    def sparse(sp_array: list):

        """
        returns a sparse array from a reduced array {sp_array}

        -------------------------------------------------------------------------
        Parameters: (Reduced array, sparsity)    -> sp_array

        -------------------------------------------------------------------------
        Examples:

        arr = ([(2, 1), (3, 2), (6, 3), (7, 1), (9, 7), (16, 2), (2, 1)], 0.67)
        sparse(arr) returns [0, 0, 1, 2 , 0, 0, 3, 1, 0, 7, 0, 0, 0, 0, 0, 0, 2, 0]

        The sparsification happens in complement to the reducing technique used in anti_sparse()

        -------------------------------------------------------------------------

        """
        array= sp_array[0]

        check= [len(i) for i in array]
        if not Metrics.equal(check):
            raise TypeError("Enter a valid sparse array")
        else:
            L, R= array[-1]
            res= [0]*L
            for i, pair in enumerate(array[:-1]):
                index= pair[0]
                value= pair[1]
                try:
                    diff= array[i+1][0]-index-1
                except:
                    diff= 0
                res.append(value)
                for i in range(diff):
                    res.append(0)

            for i in range(R):
                res.append(0)

            return res

    @staticmethod
    def integrate(f, a, b, steps=50000, precision=3):

        """
        returns the area of a function {f} from lower x bound {a} to upper x bound {b}

        """
        if abs(a-b) >= abs(steps):
            print("[DEBUG]      lower lim: {}, upper lim: {}, steps: {}".format(a, b, steps))
            print("[PROBLEM]    Steps too small for accurately calculating area.")
            print("[SOLUTION]   Please enter a higher steps or resize your domain.")
            print("[ADVICE]     The difference of the bounds must be a 1000x lesser than the steps.")
            return None

        if a == b:
            return 0
        else:
            mstart, mstop= a*steps, b*steps
            interval= (b-a)
            domain= [i/steps for i in range(mstart, mstop, interval)]
            sol_range= [f(j) for j in domain]
            w= (b-a)/steps

        return round(sum(sol_range)*w, precision)

    @staticmethod
    def get_series(f, start, num: int):

        """
        returns a list of {num} numbers which start with {start} and progress recursively using the function {f}

        -----------------------------------------------------------------------
        Parameters: function descriptor          -> f
                    initiating value             -> start
                    Number of values             -> num
        -----------------------------------------------------------------------
        Examples:

        get_series(lambda x: x**2, 3, N) returns [3, 9, 81, ... upto N values]
        get_series(lambda x: 2*x,  3, N) returns [3, 6, 12, ... upto N values]
        get_series(lambda x: x-2,  3, N) returns [3, 1, -1, ... upto N values]

        You can also use proper function definitions:

        def util(x):
            return x+2

        get_series(util, 1, N) returns N first odd numbers  [1, 3, 5, ... upto N values]
        get_series(util, 0, N) returns N first even numbers [0, 2, 4, ... upto N values]

        -------------------------------------------------------------------------

        """
        s= []
        for i in range(num):
            s.append(start)
            start= f(start)

        return s

    @staticmethod
    def equal(l: list):

        """
        returns if all the elements in a list {l} is equal numerically
        -------------------------------------------------------------------------
        Parameters: Checking List                -> l
        -------------------------------------------------------------------------
        Examples:

        equal([1])       returns True
        equal([1, 1, 1]) returns True
        equal([1, 2, 3]) returns False

        -------------------------------------------------------------------------

        """
        return len(set(l)) == 1

    @staticmethod
    def solve2D(c1: list, c2: list):

        """
        returns solutions of 2 linear equations {c1, c2} of 2 variables.
        -------------------------------------------------------------------------
        Parameters: Equation 1                   -> c1
                    Equation 2                   -> c2
        -------------------------------------------------------------------------
        Examples:

        # let equation 1 be 3x + 5y - 6= 0;
        # let equation 2 be 2x + 3y - 7= 0;

        x, y= solve2D((3, 5, -6), (2, 3, -7)) returns solutions

        -------------------------------------------------------------------------

        """

        deno  = (c1[0]*c2[1])-(c1[1]*c2[0])
        nume1 = (c1[2]*c2[1])-(c1[1]*c2[2])
        nume2 = (c1[0]*c2[2])-(c1[2]*c2[0])

        return -nume1/deno, -nume2/deno

    @staticmethod
    def strip_series(series: list):

        """
        returns list of the differences of the consecutive elements in list {series}.
        -------------------------------------------------------------------------
        Parameters: Input series                 -> series
        -------------------------------------------------------------------------
        Examples:

        series= [1, 2, 3]
        strip_series(series) returns [1, 1]

        P.S. strip_series([singular list]) returns an empty list

        """

        return [series[i]-series[i-1] for i in range(len(series))][1:]

    @staticmethod
    def get_series_depth(series: list, depth=1):

        try:
            check1= (series[1]-series[0])==(series[2]-series[1])
        except IndexError:
            raise MetricError("Low length of series")
        try:
            check2= (series[1]/series[0])==(series[2]/series[1])
        except ZeroDivisionError:
            check2= False

        if check1 or check2:
            return depth
        else:
            depth+= 1
            return Metrics.get_series_depth(Metrics.strip_series(series), depth)

    @staticmethod
    def check_return(list: list, attribute: str, val):

        """
        returns a subset of {list} where the items.attribute equals val.
        -------------------------------------------------------------------------
        Parameters: Input list                   -> list
                    Required Method              -> attribute
                    Required Value               -> val
        -------------------------------------------------------------------------
        Examples:

        list= ["apple", "ball", "cat", "dog", "ant"]

        # let's say we want a subset where the first letter is "a";
        subset= Metrics.check_return(list, "startswith", "a") returns ["apple", "ant"]

        P.S. Only attributes of list.items will work:
            startswith will work with an all string list but len will not

        """
        return [i for i in list if getattr(i, attribute)(val)]

    @staticmethod
    def rotate_array(array: list, index: int):
        """
        returns a rotated list around the specified index.
        -------------------------------------------------------------------------
        Parameters: Input list                   -> array
                    Index                        -> index
        -------------------------------------------------------------------------
        Examples:

        list= [3, 4, 7, 9, 4]

        # let's say we want a rotated array around the index 2
        {Original list:        [3, 4, 7, 9, 4]}
                                      |
        {Required return list: [9 ,4, 7, 3, 4]}

        rotated_list= Metrics.rotate_array(list, 2) returns [9, 4, 7, 3, 4]

        """
        return array[index+1:] + [array[index]] + array[:index]
