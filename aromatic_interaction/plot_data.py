import plotly.express as px
import csv


def plot_data(fname):
    with open(fname) as csvfile:
        spamreader = csv.reader(csvfile, delimiter=',')
        d=[]
        ang=[]
        info1 = []
        info2 = []
        t=0
        c=0
        for row in spamreader:
            if float(row[8])<6.0 and row[4] in ['PHE','TYR'] and row[7] in ['TRP']:
                t+=1
                #if  40.0 < float(row[9]) < 60 and 20.0 < float(row[10]) < 40.0 and 80 < float(row[12]):
                c+=1
                d.append(float(row[8]))
                ang.append(90.0-float(row[-2]))
                info1.append('{}-{}'.format(row[4],row[7]))
                info2.append('{}/{}/{}-{}-{}/{}-{}-{}'.format(row[0],row[1],row[2],row[3],row[4],row[5],row[6],row[7]))
    fig = px.scatter(x=d,y=ang,color=info1,hover_name=info2,
                     labels={'x':'Distance between ring centers','y':'Angle between mol1 ring center to mol2 ring plane'})
    fig.add_vrect(x0=4.7, x1=6.7)
    fig.add_hline(y=5)
    fig.add_hline(y=15)
    fig.show()
    fig.write_html('Distance_vs_angle3.html')
    print ('Total {}'.format(t))
    print ('Cross linked {}'.format(c))


if __name__ == "__main__":
    plot_data('data/pdb_aro.csv')