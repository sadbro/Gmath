import numpy as np
from shapely.geometry import LineString as LS
import matplotlib.pyplot as plt
from metrics import *

__all__ =   {
                "pi":Metrics.pi,
                "e":Metrics.e,
                "sin":Metrics.sin,
                "cos":Metrics.cos,
                "log":Metrics.log,
                "sqrt":Metrics.sqrt,
                "fact":Metrics.fact
            }

class Gfunc:

    def __init__(self, expression_string, var="x", start=0, end=10, steps=1000, optx=None, inclusive=True, local_constants={}):

        self.expression= expression_string
        self.input_var= var
        if optx is None:
            self.title= "f({})= {}".format(self.input_var, self.expression)
        else:
            self.title= "f({})= {}".format(self.input_var, optx)
        self.start_x= start
        self.end_x= end
        self.steps_x= steps
        self.width_x= (self.end_x-self.start_x)/self.steps_x
        self.locales= local_constants
        if inclusive:
            self.domain= list(np.linspace(self.start_x, self.end_x, self.steps_x))
        else:
            self.domain= list(np.linspace(self.start_x+self.width_x, self.end_x, self.steps_x))
        self.range= self.function(self.domain)
        self.xy= self.domain, self.range
        self.yx= (self.xy)[::-1]

    def __str__(self):

        object= "\nOBJECT: geom.Gfunc"
        func= "\nFUNCTION: {}".format(self.title)
        in_var= "\nINPUT VARIABLE: {}".format(self.input_var)
        domain= "\nDOMAIN: {}".format([round(self.start_x, 4), round(self.end_x, 4)])
        fr= np.array(self.range)
        r= list(fr[~np.isnan(fr)])
        range= "\nRANGE: {}\n".format([round(min(r), 4), round(max(r), 4)])

        return object + func + in_var + domain + range

    def __add__(self, other):

        if self.input_var != other.input_var:
            raise TypeError("Objects are not of same variable")

        else:
            return Gfunc(self.expression+"+"+other.expression, var=self.input_var, start=self.start_x, end=self.end_x, steps=self.steps_x)

    def __sub__(self, other):

        if self.input_var != other.input_var:
            raise TypeError("Objects are not of same variable")

        else:
            return Gfunc(self.expression+"-"+other.expression, var=self.input_var, start=self.start_x, end=self.end_x, steps=self.steps_x)

    def __mul__(self, other):

        if self.input_var != other.input_var:
            raise TypeError("Objects are not of same variable")

        else:
            return Gfunc(self.expression+"*"+other.expression, var=self.input_var, start=self.start_x, end=self.end_x, steps=self.steps_x)

    def __truediv__(self, other):

        if self.input_var != other.input_var:
            raise TypeError("Objects are not of same variable")

        else:
            return Gfunc(self.expression+"/"+other.expression, var=self.input_var, start=self.start_x, end=self.end_x, steps=self.steps_x)

    def compare(self, other):

        plt.plot(*self.xy, label=self.expression)
        plt.plot(*other.xy, label=other.expression)

        for x, y in zip(*self.intersect(other)):
            plt.plot(x, y, 'ro')

        plt.title("[{}] ~ [{}]".format(self.expression, other.expression))
        plt.legend(shadow=True)
        plt.grid(True)
        plt.show()

        return self.intersect(other)

    def area_intersect(self, other, start, end):

        self_area= self.find_area(start, end)
        other_area= other.find_area(start, end)

        return abs(self_area-other_area)

    def intersect(self, other, domain_check=False):

        if domain_check:
            if self.domain != other.domain:
                raise ValueError("Conflicting Domains found!!")

        self_ls= LS(np.column_stack((self.xy)))
        other_ls= LS(np.column_stack((other.xy)))

        x_coord= []
        y_coord= []
        ittn= self_ls.intersection(other_ls)

        try:
            try:
                for p in ittn:
                    x_coord.append(p.x)
                    y_coord.append(p.y)

            except TypeError:
                x_coord.append(ittn.x)
                y_coord.append(ittn.y)

        except AttributeError:
            pass

        return x_coord, y_coord

    def function(self, input_value, r=5):

        if type(input_value) == list:
            output= []
            for i in input_value:
                exec("{}= {}".format(self.input_var, i))
                self.locales[self.input_var]= i
                try:
                    output.append(round(eval(self.expression, self.locales), r))
                except ValueError:
                    output.append(np.nan)
                except TypeError:
                    output.append(np.nan)
                except ZeroDivisionError:
                    output.append(np.inf)
                except OverflowError:
                    output.append(np.nan)
        else:
            exec("{}= {}".format(self.input_var, input_value))
            self.locales[self.input_var]= input_value
            try:
                output= round(eval(self.expression, self.locales), r)
            except ValueError:
                output= np.nan
            except ZeroDivisionError:
                output= np.inf

        return output

    @staticmethod
    def show_multiple(*Gfunc_objects):

        for function in Gfunc_objects:
            plt.plot(function.domain, function.range, label=function.title)
        plt.grid(True)
        plt.legend(loc="best")
        plt.show()

    def find_area(self, a, b, steps=50000, precision=3):

        return Metrics.integrate(self.function, a, b, steps, precision)

    def embed(self, other):

        return Gfunc(self.expression.replace(self.input_var, other.expression), var=other.input_var)

    def show(self, is_centered=False):

        print(self)

        x_output= self.domain
        y_output= self.range

        if is_centered:

            fig, axes= plt.subplots()

            axes.spines['left'].set_position('center')
            axes.spines['bottom'].set_position('center')

            axes.spines['top'].set_color('none')
            axes.spines['right'].set_color('none')

            axes.xaxis.set_ticks_position('bottom')
            axes.yaxis.set_ticks_position('left')

        plt.title(self.title)
        plt.plot(x_output, y_output)
        plt.grid(True)
        plt.show()

    def find(self, y=0, r=4):

        xy= self.domain, [y]*self.steps_x
        self_ls= LS(np.column_stack((self.xy)))
        other_ls= LS(np.column_stack((xy)))

        x_coord= []
        y_coord= []
        try:
            ittn= self_ls.intersection(other_ls)
        except:
            raise ValueError()

        try:
            try:
                for p in ittn:
                    x_coord.append(p.x)
                    y_coord.append(p.y)

            except TypeError:
                x_coord.append(ittn.x)
                y_coord.append(ittn.y)

            except AttributeError:
                x_coord= [i.xy[0][0] for i in ittn]
                y_coord= [i.xy[1][0] for i in ittn]

        except:
            raise AttributeError("NO INTERSECTION FOUND AT Y = {} for {}".format(y, self.title))

        plt.plot(*self.xy, label=self.title)
        plt.plot(*xy, label="{}".format(y))

        for x, y in zip(x_coord, y_coord):
            plt.plot(x, y, 'ro')

        plt.title("{} = {}".format(self.title, y))
        plt.legend(shadow=True)
        plt.grid(True)
        plt.show()

        res= []
        for x, y in zip(x_coord, y_coord):
            x= round(x, r)
            res.append((x, y))

        return res

    def t_find(self, y=0, r=4):

        for j, i in enumerate(self.range):
            if i == round(y, 4):
                return round(self.domain[j]), y

    def differentiate(self, only_return=False):

        dys= []
        for index, i in enumerate(self.domain):
            try:
                dy= (self.range[index+1]-self.range[index])/self.width_x
                dys.append(dy)
                #print("[{}, {}]".format(i, self.range[index]))
            except:
                pass

        dys.append(dys[-1])

        if only_return:
            return dys
        else:
            plt.plot(*self.xy, label=self.title)
            plt.plot(self.domain, dys, label="derivative")
            plt.legend(loc="best")
            plt.grid()
            plt.show()

    def der_find(self, y=0, show_og=True, only_return=False):

        xy= self.domain, [y]*self.steps_x
        self_ls= LS(np.column_stack((self.domain, self.differentiate(True))))
        other_ls= LS(np.column_stack((xy)))

        x_coord= []
        y_coord= []
        try:
            ittn= self_ls.intersection(other_ls)
        except:
            raise ValueError()

        try:
            try:
                for p in ittn:
                    x_coord.append(p.x)
                    y_coord.append(p.y)

            except TypeError:
                x_coord.append(ittn.x)
                y_coord.append(ittn.y)

            except AttributeError:
                x_coord= [i.xy[0][0] for i in ittn]
                y_coord= [i.xy[1][0] for i in ittn]

        except:
            raise AttributeError("NO INTERSECTION FOUND AT Y = {} for f'({})".format(y, self.input_var))

        res= []
        for x, y in zip(x_coord, y_coord):
            rx= round(x, 3)
            res.append((rx, y))

        if only_return:
            return res

        else:
            if show_og:
                plt.plot(*self.xy, label=self.title)

            plt.plot(self.domain, self.differentiate(True), label="f'({})".format(self.input_var))

            for x, y in zip(x_coord, y_coord):
                plt.plot(x, y, 'ro')

            plt.title("{} = {}".format("f'({})".format(self.input_var), y))
            plt.legend(shadow=True)
            plt.grid(True)
            plt.show()

    @staticmethod
    def parametric(xg1, yg2, aspect_eq=True):

        if len(xg1.range) != len(yg2.range):
            raise MetricError("Unequal dimensions for parametrization: {} ~ {}".format(len(xg1.range), len(yg2.range)))
        plt.plot(xg1.range, yg2.range)
        plt.title("y={}, x={}".format(yg2.expression, xg1.expression))
        if aspect_eq:
            plt.gca().set_aspect("equal", adjustable="box")
        plt.grid(True)
        plt.show()

g1 = Gfunc("sin(x)", local_constants=__all__)

g1.find(0)
