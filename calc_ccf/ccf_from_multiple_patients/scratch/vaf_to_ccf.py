import numpy
def vaf_to_ccf(vaf,pur,prev,minor_cn,major_cn):
    cn = major_cn + minor_cn
    alpha = (cn*prev + 2*(1-prev))*pur + 2*(1-pur)
    #no duplications
    if (minor_cn <= 1 & major_cn <= 1):
        return alpha*vaf/pur
    else:
        #One duplication, no LOH
        if (minor_cn == 1):
            if (vaf >= (major_cn*prev + 1)*pur/2*alpha):
                return (alpha*vaf/pur)-(prev*(major_cn-1))
            else:
                return alpha*vaf/pur
        #One duplication with LOH
        else:
            if(minor_cn==0):
                if(vaf >= (1+(major_cn-1)*prev)*pur/(2*alpha)):
                   return alpha*vaf/pur - (major_cn-1)*prev
                else:
                    return alpha*vaf/pur
            #two duplications
            else:
                if(vaf <= (1+(minor_cn-1)*prev)*pur/(2*alpha)):
                   return alpha*vaf/pur
                else:
                   if(vaf >= (1+(major_cn+minor_cn-1)*prev)*pur/(2*alpha)):
                      return alpha*vaf/pur - (major_cn-1)*prev
                   else:
                      return alpha*vaf/pur - (minor_cn-1)*prev

