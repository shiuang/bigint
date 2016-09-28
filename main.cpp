#include <stdio.h>
#include <stdlib.h>

//#define FILEIO
#ifdef FILEIO
#include <time.h>
#endif

const int MAXN = 1<<17;
struct bigint
{
	int data[MAXN+10];
	int length;
	int sign;
}first,second,ans,ans1;
char op;

void copy(const int *x, int *y, int l);

char input[MAXN];

int read(bigint &input_num);
int read_char(char &op);
int read_enter();
void write(const bigint &a, const char end);

int cmp(const bigint &a, const bigint &b, int st = 0);
void plus(bigint &a, bigint &b, bigint &c);
void minus_t(const bigint &a, const bigint &b, bigint &c, int st = 0);
void minus(bigint &a, bigint &b, bigint &c);


typedef long long ll;
const int P = (479<<21)+1;
int g;
int w[2][MAXN];
int tmp[MAXN];

void init_NTT();
int power(ll a, int b, int mod);
void swap(int &a, int &b);
void NTT(int *x, int l, int on);
void mul(const bigint &a, const bigint &b, bigint &c);

const int sum = 5000;
bigint tmp1,two;

void div(const bigint &a, const bigint &b, bigint &c, bigint &d);

void copy(const int *x, int *y, int l)
{
	int i;
	for (i = 0; i<l; i++)
		y[i] = x[i];
}

//read&write

int read(bigint &input_num)
{
	char c,flag=0;
	int i = 0;
	input_num.sign = 0;
	while (~(c = getchar()))
	{
		if (c=='\n') break;
		if (c>='0' && c<='9')
		{
			if (flag==4) return 1;
			if (c!='0' || flag==3)
			{
				if (i==100000) return 1;
				flag = 3;
				input[i++] = c - '0';
			}
			else flag = 2;
		}
		else if (c=='+' || c=='-')
		{
			if (flag!=0) return 1;
			flag = 1;
			if (c=='-') input_num.sign = 1;
		}
		else if (c==' ')
		{
			if (flag>=2)
				flag = 4;
			else if (flag!=0)
				return 1;
		}
		else return 1;
	}
	if (c==-1)
	{
		if (flag==2 || flag==1)
			return 2;
		else 
			return 3;
	}
	if (flag<2) return 4;

	int j=0,l=0;
	for (i--; i>=0; l++, i-=2)
	{
		input_num.data[l] = input[i];
		if (i-1>=0)
			input_num.data[l] += input[i-1]*10;
	}
	input_num.length = l;

	return 0;
}

int read_char(char &op)
{
	int flag = 0;
	char c;
	while (~(c = getchar()))
	{
		if (c=='\n') break;
		switch (c)
		{
		case '+':case '-':case '*':case '/':
			if (flag!=0) return 1;
			flag = 1;
			op = c;
			break;
		case ' ':
			if (flag==1) flag = 2;
			break;
		default:
			return 1;
		}
	}
	if (c==-1)
	{
		if (flag==0)
			return 2;
		else return 3;
	}
	if (flag==0) return 4;
	return 0;
}

int read_enter()
{
	char c;
	while (~(c = getchar())) if (c=='\n') return 1;
	return 2;
}

void write(const bigint &a, const char end)
{
	if (a.length==0)
	{
		printf("0%c",end);
		return;
	}

	int i;
	if (a.sign)
		putchar('-');
	for (i = a.length-1; i>=0; i--)
	{
		if (a.data[i]<=9 && i!=a.length-1) putchar('0');
		printf("%d",a.data[i]);
	}
	putchar(end);
}

//plus&minus

void plus(bigint &a, bigint &b, bigint &c)
{
	if (a.length==0 && b.length==0)
	{
		c.length = 0;
		return;
	}

	if (a.sign!=b.sign)
	{
		if (a.sign)
		{
			a.sign = 0;
			minus(b,a,c);
		}
		else
		{
			b.sign = 0;
			minus(a,b,c);
		}
		return;
	}

	c.sign = a.sign;
	c.length = a.length;
	if (b.length>c.length) c.length = b.length;

	int i;
	c.data[0] = 0;
	for (i = a.length; i<c.length; i++)
	{
		a.data[i] = 0;
	}
	for (i = b.length; i<c.length; i++)
	{
		b.data[i] = 0;
	}
	for (i = 0; i<c.length; i++)
	{
		c.data[i] += a.data[i]+b.data[i];
		if (c.data[i]<100)
			c.data[i+1] = 0;
		else
		{
			c.data[i+1] = 1;
			c.data[i] -= 100;
		}
	}
	if (c.data[c.length]) c.length++;
}

int cmp(const bigint &a, const bigint &b, int st)
{
	if (a.length-st>b.length) return 1;
	if (a.length-st<b.length) return -1;

	int i;
	for (i = b.length-1; i>=0; i--)
	{
		if (a.data[i+st]>b.data[i]) return 1;
		if (a.data[i+st]<b.data[i]) return -1;
	}
	return 0;
}

void minus(bigint &a, bigint &b, bigint &c)
{
	if (a.length==0 && b.length==0)
	{
		c.length = 0;
		return;
	}

	if (a.sign!=b.sign)
	{
		b.sign = a.sign;
		plus(a,b,c);
		return;
	}

	int i;
	for (i = a.length; i<b.length; i++)
	{
		a.data[i] = 0;
	}
	for (i = b.length; i<a.length; i++)
	{
		b.data[i] = 0;
	}

	int sign = cmp(a,b);
	if (sign==0)
	{
		c.length = 0;
		return;
	}
	if (sign>0)
	{
		c.sign = a.sign;
		minus_t(a,b,c);
	}
	else
	{
		c.sign = 1-a.sign;
		minus_t(b,a,c);
	}
}

void minus_t(const bigint &a, const bigint &b, bigint &c, int st)
{
	int i;
	c.length = a.length;
	c.data[st] = 0;
	for (i = st; i<c.length; i++)
	{
		c.data[i] += a.data[i] - b.data[i-st];
		if (c.data[i]<0)
		{
			c.data[i+1] = -1;
			c.data[i] += 100;
		}
		else
			c.data[i+1] = 0;
	}
	while (c.data[c.length-1]==0) c.length--;
}

//mul
//NTT算法
//http://www.cnblogs.com/Trinkle/p/3856040.html

void init_NTT()
{
	int i;
	g = power(3,(P-1)/MAXN,P);
	w[0][0] = w[1][0] = 1;
	for (i = 1; i<MAXN; i++)
		w[0][i] = w[1][MAXN-i] = (ll)w[0][i-1]*g%P;
}

int power(ll a, int b, int mod)
{
	ll ans = 1;
	while (b)
	{
		if (b&1) ans = ans*a%mod;
		a = a*a%mod;
		b >>= 1;
	}
	return ans;
}

void swap(int &a, int &b)
{
	int tmp;
	tmp = a;
	a = b;
	b = tmp;
}

void NTT(int *x, int l, int on)
{
	int i,j,k;
	for (i = j = 0; i<l; i++)
	{
		if (i>j) swap(x[i],x[j]);
		for (k = l>>1; j&k; k>>=1) j ^= k;
		j |= k;
	}

	int tmp;
	for (i = 1; i<=l>>1; i<<=1)
		for (j = 0; j<l; j+=i<<1)
			for (k = 0; k<i; k++)
			{
				tmp = (ll)x[j+k+i]*w[on][MAXN/(i<<1)*k]%P;
				x[j+k+i] = (x[j+k]-tmp+P)%P;
				x[j+k] = (x[j+k]+tmp)%P;
			}
}

void mul(const bigint &a, const bigint &b, bigint &c)
{
	if (a.length==0 || b.length==0)
	{
		c.length = 0;
		return;
	}
	c.sign = a.sign^b.sign;
	c.length = a.length + b.length - 1;
	
	int i;
	int length2;
	for (length2 = 1; length2<c.length; length2<<=1);
	copy(a.data, tmp, a.length);
	copy(b.data, c.data, b.length);
	for (i = a.length; i<length2; i++)
		tmp[i] = 0;
	for (i = b.length; i<=length2; i++)
		c.data[i] = 0;

	NTT(tmp, length2, 0);
	NTT(c.data, length2, 0);
	
	for (i = 0; i<length2; i++)
		c.data[i] = (ll)tmp[i]*c.data[i]%P;

	NTT(c.data, length2, 1);

	int inv = power(length2,P-2,P);
	for (i = 0; i<c.length; i++)
		c.data[i] = (ll)c.data[i]*inv%P;
	for (i = 0; i<c.length; i++)
	{
		c.data[i+1] += c.data[i]/100;
		c.data[i] %= 100;
	}
	if (c.data[c.length]) c.length++;
}

//div
//牛顿迭代法x(n+1)=x(n)*(2-b*x(n))计算1/b
//通过计算解的前几位来估计初值x(0)

void div(const bigint &a, const bigint &b, bigint &c, bigint &d)
{
	if (a.length==0)
	{
		c.length = d.length = 0;
		return;
	}
	c.sign = d.sign = 0;

	int sign = cmp(a,b);
	int i;
	if (sign<0)
	{
		c.length = 0;
		d.length = a.length;
		copy(a.data, d.data, a.length);
		return;
	}
	if (sign==0)
	{
		c.length = c.data[0] = 1;
		d.length = 0;
		return;
	}

	int iter=sum/b.length+1;
	tmp1.sign = 0;
	tmp1.length = a.length + 2;
	if (iter+b.length<tmp1.length)
		tmp1.length = iter+b.length;
	for (i = 0; i<tmp1.length-1; i++)
		tmp1.data[i] = 0;
	tmp1.data[i] = 1;
	c.length = tmp1.length - b.length;
	for (i = c.length-1; i>=0; i--)
	{
		c.data[i] = 0;
		while (cmp(tmp1,b,i)>=0)
		{
			d.length = tmp1.length;
			copy(tmp1.data, d.data, tmp1.length);
			minus_t(d, b, tmp1, i);
			c.data[i]++;
		}
	}
	
	tmp1.length = a.length + 2 - b.length;
	iter = tmp1.length - c.length;
	for (i = 0; i<iter; i++) tmp1.data[i] = 0;
	for (; i<tmp1.length; i++) tmp1.data[i] = c.data[i-iter];
	if (tmp1.data[i-1]>=100)
	{
		tmp1.data[i-1] -= 100;
		tmp1.data[tmp1.length++] = 1;
	}

	two.sign = 0;
	two.length = a.length+2;
	for (i = 0; i<two.length-1; i++)
		two.data[i] = 0;
	two.data[i] = 2;
	c.length = 0;
#ifdef FILEIO
		iter = 0;
#endif
	while (cmp(tmp1,c)!=0)
	{
#ifdef FILEIO
		iter++;
#endif

		c.length = tmp1.length;
		copy(tmp1.data, c.data, tmp1.length);

		mul(b,c,tmp1);
		for (i = tmp1.length; i<two.length; i++) tmp1.data[i] = 0;
		minus_t(two,tmp1,d);
		mul(c,d,tmp1);

		for (i = 0; i+a.length+1<tmp1.length; i++)
			tmp1.data[i] = tmp1.data[i+a.length+1];
		tmp1.length -= a.length+1;
	}
#ifdef FILEIO
	printf("%d\n",iter);
#endif

	mul(tmp1,a,c);
	for (i = 0; i+a.length+1<c.length; i++)
		c.data[i] = c.data[i+a.length+1];
	c.length -= a.length+1;
	mul(b,c,tmp1);
	tmp1.data[tmp1.length] = 0;
	minus_t(a,tmp1,d);

	sign = cmp(d,b);
	if (sign>=0)
	{
		for (i = 0; i<c.length && c.data[i]==99; i++)
			c.data[i]=0;
		if (i==c.length)
			c.data[c.length++] = 0;
		c.data[i]++;
		if (sign==0)
		{
			d.length = 0;
		}
		else
		{
			tmp1.length = d.length;
			copy(d.data, tmp1.data, d.length);
			minus_t(tmp1,b,d);
		}
	}
}

int main()
{
#ifdef FILEIO
	freopen("input.txt","r",stdin);
	freopen("output.txt","w",stdout);
	clock_t start,end;
#endif
	
	init_NTT();
	int status;
	while (1)
	{
		status = read(first);
		if (status==1) status = read_enter();
		if (status==1 || status==2 || status==4)
		{
			printf("Error input\n");
			if (status==1 || status==4) status = read_enter();
			if (status==1) status = read_enter();
		}
		if (status==2 || status==3) break;
		if (status==1) continue;

		status = read(second);
		if (status==1) status = read_enter();
		if (status==1 || status==2 || status==4)
		{
			printf("Error input\n");
			if (status==1 || status==4) status = read_enter();
		}
		if (status==2 || status==3) break;
		if (status==1) continue;

		status = read_char(op);
		if (status==1) status = read_enter();
		if (status==1 || status==2 || status==4)
		{
			printf("Error input\n");
		}
		if (status==2) break;
		if (status==1 || status==4) continue;

		switch (op)
		{
		case '+':
			plus(first,second,ans);
			write(ans,'\n');
			break;
		case '-':
			minus(first,second,ans);
			write(ans,'\n');
			break;
		case '*':
#ifdef FILEIO
			start=clock();
#endif
			mul(first,second,ans);
#ifdef FILEIO
			end=clock();
			printf("%d\n",end-start);
#endif
			write(ans,'\n');
			break;
		case '/':
			if ((first.sign && first.length!=0) || second.sign || second.length==0)
			{
				status = 1;
				break;
			}
			first.sign = 0;
			second.data[second.length] = 0;
#ifdef FILEIO
			start=clock();
#endif
			div(first,second,ans,ans1);
#ifdef FILEIO
			end=clock();
			printf("%d\n",end-start);
#endif
			write(ans,' ');
			write(ans1,'\n');
			break;
		}

		if (status==1) printf("Error input\n");
		if (status==3) break;
	}
	fclose(stdin);
	fclose(stdout);
	system("pause");
	return 0;
}