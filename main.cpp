#include <iostream>
#include <cmath>
#include <ctime>
#include<cstdlib>

using namespace std;

// 定义数组：colony[][]:种群数组   fitness[]：种群适应值  Distance[]种群距离  citydistance[]：各城市距离  best_path[]:最佳种群
// 定义变量：city_num：基因个数    pop_size:种群个数   improve_count:改良次数   itter_time:进化代数
// 相关概率： pm:变异概率    pc:交叉概率
const int city_num = 10;
const int pop_size = 20;
const int improve_time = 5;
const int itter_time = 20;
const int start = 0;
const float pm = 0.4;
const float pc = 0.8;
const int N = 10000;

double city_distance[city_num][city_num];
double city_condition[city_num][2];

typedef struct {
    int colony[pop_size][city_num + 1];
    double fitness[pop_size];
    double Distance[pop_size];
    double bestfitness;
    double bestdistance;
    int bestroot;
    int bestrooting[city_num + 1];
} TSP, *P;              //进行TSP分析的基本单位


void get_condition();
//获得每个城市的坐标
void get_distance();
//获得各个城市间的距离
int check(TSP &city, int pop_test, int k, int city_test);
//检查产生的基因是否在染色体内已经出现
void initial(TSP &city);
//产生初代种群
void improve(TSP &city, int pop_num, double before);
//对初代种群进行事先优化
void copy(int target[], int source[]);
//复制数组以便备份数据
double get_one_distance(int a[]);
//获得路径距离
void disp(TSP &city);
//输出最优路径和距离
void calculate(TSP &city);
//更新每个种群的各项数据（最佳路径，每条路径距离，适应值)
void natural_select(TSP &city);
//自然选择
void mutation(TSP &city);
//基因变异
void disp_result(TSP &city);
//输出所有染色体
void crossover(TSP &city);
//基因重组
int check_cross(int *p, int k, int num);
//检查是否出现（同check（））

void get_condition() {
    int i = 0;
    double x, y;
    for (i = 0; i < city_num; i++) {
        cout << "input no." << i << ".x"<<endl;
        cin >> x;
        cout << "input no." << i << ".y"<<endl;
        cin >> y;
        city_condition[i][0] = x;
        city_condition[i][1] = y;           //保存各城市数据
    }
   // for (i = 0; i < city_num; i++) {
        //cout << "city " << i + 1 << endl;
        //cout << "(" << city_condition[i][0] << "," << city_condition[i][1] << ")" << endl;
    //}                                       //输出各城市坐标
}

void get_distance() {
    int i, j;
    double x, y;
    for (i = 0; i < city_num; i++) {
        city_distance[i][i] = 0;
        for (j = i + 1; j < city_num; j++) {
            x = city_condition[i][0] - city_condition[j][0];
            y = city_condition[i][1] - city_condition[j][1];
            city_distance[i][j] = sqrt(x * x + y * y);          //计算距离
        }
    }
    for (i = 0; i < city_num; i++) {
        for (j = 0; j < i; j++) {
            city_distance[i][j] = city_distance[j][i];
        }
    }
    //cout << "city_distance:" << endl;
    //for (i = 0; i < city_num; i++) {
        //for (j = i + 1; j < city_num; j++) {
            //cout << "no." << i << " to no." << j << " is " << city_distance[i][j] << endl;
        //}
    //}                                                             //输出各城市距离
}

int check(TSP &city, int pop_test, int k, int city_test) {
    int i;
    for (i = 0; i < city_test; i++) {
        if (city.colony[pop_test][i] == k)
            return false;
    }
    return true;
}

void initial(TSP &city) {
    int i, j, r = 0, flag;
    for (i = 0; i < pop_size; i++) {
        city.colony[i][0] = start;
        city.colony[i][city_num] = start;
        city.Distance[i] = 0;
        city.fitness[i] = 0;
    }
    city.bestdistance = 0;
    city.bestfitness = 0;
    city.bestroot = 0;
    for (i = 0; i < city_num + 1; i++)
        city.bestrooting[i] = 0;                        //初始化
    //cout << "ok" << endl;
    for (i = 0; i < pop_size; i++) {
        for (j = 1; j < city_num; j++) {
            flag = 0;
            while (!flag) {
                r = rand() % (city_num);                //随机产生基因
                flag = check(city, i, r, j);         //检查是否重复
            }
            city.colony[i][j] = r;                      //加入初始路径
        }
    }
    //cout << "before improve:" << endl;
    //disp(city);
    double before;
    for (i = 0; i < pop_size; i++) {
        before = get_one_distance(city.colony[i]);
        improve(city, i, before);                   //改进随机生成的路径
    }
    //cout << "after improve:" << endl;
    disp(city);
}

void disp(TSP &city) {
    cout << "======the root=====" << endl;
    int i, j;
    for (i = 0; i < pop_size; i++) {
        cout << "Rooting " << i << endl;
        for (j = 0; j < city_num + 1; j++) {
            cout << city.colony[i][j] << "->";
        }
        cout << endl;
    }                                               //输出生成的种群
    cout << "=====end=====" << endl;
}

void improve(TSP &city, int pop_num, double before) {
    int a[city_num + 1], b[city_num + 1];
    int *p = city.colony[pop_num];
    int i = 0, change_1, change_2, t;
    double after;
    copy(a, p);                                 //备份
    while (i < improve_time) {
        change_1 = rand() % (city_num);
        change_2 = rand() % (city_num);
        if (change_1 && change_2) {             //随机选择内部基因
            if (change_1 != change_2) {
                copy(b, a);
                t = b[change_1];
                b[change_1] = b[change_2];
                b[change_2] = t;                //交换两处基因
                after = get_one_distance(b);    //获得交换后距离
                if (after < before) {
                    copy(a, b);
                    before = after;             //如果路径变小，则替换原路径
                }
            }
        }
        i++;
    }
    copy(p, a);                                 //如果路径没变更优，则不改变
}

void copy(int target[], int source[]) {
    for (int i = 0; i < city_num + 1; i++)
        target[i] = source[i];
}

double get_one_distance(int a[]) {
    int i;
    double distance = 0;
    for (i = 0; i < city_num; i++) {
        distance += city_distance[a[i]][a[i + 1]];      //获得指定路径的距离
    }
    return distance;
}

void calculate(TSP &city) {
    int i;
    for (i = 0; i < pop_size; i++) {
        city.Distance[i] = get_one_distance(city.colony[i]);
        city.fitness[i] = 1.0 / city.Distance[i];
    }                                                   //计算种群距离和适应值
    city.bestfitness = city.fitness[0];
    city.bestdistance = city.Distance[0];
    for (i = 0; i < pop_size; i++) {
        if (city.fitness[i] > city.bestfitness)
            city.bestfitness = city.fitness[i];           //获得best_fitness
        if (city.Distance[i] < city.bestdistance)
            city.bestdistance = city.Distance[i];         //获得best_distance
    }
    for (i = 0; i < pop_size; i++) {
        if (city.Distance[i] == city.bestdistance) {
            city.bestroot = i;
            copy(city.bestrooting, city.colony[i]);
        }                                               //获得bestrooting[city_size]
    }
}

void natural_select(TSP &city) {
    double sum_fitness = 0, rate;
    double p[pop_size];
    double cumulative_P[pop_size];
    int temp[pop_size][city_num + 1] = {0};
    int i, j;
    copy(temp[0], city.colony[city.bestroot]);       //保存最佳路径
    for (i = 0; i < pop_size; i++)
        sum_fitness += city.fitness[i];                     //sum_fitness（总适应值）
    for (i = 0; i < pop_size; i++)
        p[i] = city.fitness[i] / sum_fitness;               //p(概率)
    cumulative_P[0] = p[0];
    for (i = 0; i < pop_size - 1; i++)
        cumulative_P[i + 1] = cumulative_P[i] + p[i + 1];     //cumulative_P（累计概率）
    for (i = 1; i < pop_size; i++) {
        rate = rand() * 1.0 / RAND_MAX;                     //随机概率
        //cout << "select_rate:" << rate << endl;
        for (j = 0; j < pop_size; j++) {
            if (cumulative_P[j] >= rate)
                break;                                      //选中
        }
        //cout << "choose No." << j << endl;
        copy(temp[i], city.colony[j]);                      //选中则进入temp数组存储
    }
    for (i = 0; i < pop_size; i++)
        copy(city.colony[i], temp[i]);                      //替换原种群
    //disp(city);
}

void mutation(TSP &city) {
    int i, left, right, j, t, r;
    double rate_mutation;
    int temp[city_num + 1];
    double after, before;
    for (i = 0; i < pop_size; i++) {
        if (i != city.bestroot) {                                  //保留最佳路径
            rate_mutation = rand() * 1.0 / RAND_MAX;               //根据pm判断是否基因重组
            //cout << "rate_mutation:" << rate_mutation << endl;
            if (rate_mutation > pm) {
                left = rand() % (city_num);
                right = rand() % (city_num);                       //随机选择两处基因
                //cout << "left:" << left << endl;
                //cout << "right:" << right << endl;
                if (left != 0 && right != 0) {
                    if (abs(left - right) <= 1) {              //避免无效交换
                        //cout << "false,next_time;" << endl;
                        continue;                                  //保证不选择出发城市
                    } else if (left > right) {
                        t = left;
                        left = right;
                        right = t;
                    }                                              //保证left<right
                    before = get_one_distance(city.colony[i]);
                    //cout << "before:" << before << endl;
                    copy(temp, city.colony[i]);                  //备份
                    t = (left + right) / 2;
                    for (j = left; j < t; j++) {
                        r = temp[j];
                        temp[j] = temp[2 * t - j];
                        temp[2 * t - j] = r;
                    }                                           //将选中基因片段首尾两两交换
                    after = get_one_distance(temp);
                    //cout << "after:" << after << endl;
                    if (before > after) {
                        copy(city.colony[i], temp);
                        //cout << "change;" << endl;          //如果变异后更优，则替换原染色体
                    } //else
                        //cout << "not change;" << endl;
                }
            } //else
                //cout << "next_time;" << endl;
        }
    }
    //disp(city);
}

void disp_result(TSP &city) {
    int i;
    //for (i = 0; i < pop_size; i++)
        //cout << "No." << i << ":" << city.Distance[i] << endl;
    //cout << "the best root is:" << endl;
    for (i = 0; i < city_num + 1; i++) {
        cout << city.bestrooting[i] << "->";
    }
    cout << endl;
    cout << "the minest distance is:" << city.bestdistance << endl;
    //cout << "the largest fitness is:" << city.bestfitness << endl;
}

void crossover(TSP &city) {
    int i, j, mom, dad, flag = 0, k;
    int temp1[city_num + 1] = {0};
    double cross_rate, before, after;
    for (i = 0; i < pop_size; i++) {
        //cout << endl;
        //cout << "round " << i << endl;
        if (i != city.bestroot) {                           //保留最佳路径
            cross_rate = rand() * 1.0 / RAND_MAX;
            //cout << "cross_rate:" << cross_rate << endl;
            if (cross_rate < pc) {                          //通过pc控制是否重组
                mom = rand() % (pop_size);
                dad = rand() % (pop_size);                  //随机选择父母染色体
                //cout << "mom:" << mom << endl;
                //cout << "dad:" << dad << endl;
                if (mom == dad)
                    continue;                               //避免无效重组
                copy(temp1, city.colony[i]);                //备份
                //cout << "mom:" << endl;
                //for (j = 0; j < city_num + 1; j++)
                    //cout << city.colony[mom][j] << "->";
                //cout << endl;
                //cout << "dad:" << endl;
                //for (j = 0; j < city_num + 1; j++)
                    //cout << city.colony[dad][j] << "->";
                //cout << endl;
                for (j = 1; j < city_num; j++) {
                    if (j % 2 == 0)
                        temp1[j] = city.colony[mom][j];
                    else
                        temp1[j] = city.colony[dad][j];
                }                                           //将父母段交替赋给子段
                //cout << "before cross:" << endl;
                //for (j = 0; j < city_num + 1; j++)
                    //cout << temp1[j] << "->";
                //cout << endl;
                for (j = 0; j < city_num; j++) {
                    flag = check_cross(temp1, temp1[j], j);        //检查是否重复
                    k = 0;
                    while (flag) {
                        temp1[j] = city.colony[mom][k];
                        flag = check_cross(temp1, temp1[j], j);
                        k++;
                    }                                           //若重复，则到父段中寻找，直至不再重复
                }
                before = get_one_distance(city.colony[i]);
                after = get_one_distance(temp1);
                //cout << "after cross:" << endl;
                //for (j = 0; j < city_num + 1; j++)
                    //cout << temp1[j] << "->";
                //cout << endl;
                if (before > after){                            //若重组后更优，则替换
                    copy(city.colony[i], temp1);
                    //cout<<"cross!"<<endl;
                }
                 //else
                    //cout << "not cross!" << endl;
            }
        }
    }
    //disp(city);
}

int check_cross(int *p, int k, int num) {
    for (int i = 0; i < num; i++) {
        if (*(p + i) == k)
            return 1;
    }
    return 0;
}

int main() {
    TSP city;
    int i;
    srand((int) time(0));
    get_condition();                                    //获得坐标
    get_distance();                                     //获得距离
    initial(city);                                  //初始化种群并优化
    calculate(city);                                //计算相关数据
    //disp_result(city);                              //展示初代种群
    for (i = 0; i < itter_time; i++) {
        cout << endl;
        cout << "====next:" << i << "=====" << endl;
        natural_select(city);                       //自然选择
        mutation(city);                             //基因变异
        crossover(city);                            //基因重组
        calculate(city);
        disp_result(city);
    }
    //disp_result(city);
}




