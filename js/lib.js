var getParams = function($form) {
  params = {}
  $.each($form.serializeArray(), function() {
    params[this.name] = this.value;
  });

  params['gamma'] = params['gamma'].split(' ');
  params['epsilon'] = params['epsilon'].split(' ');
  params['lambda'] = params['lambda'].split(' ');
  params['phi'] = params['phi'].split(' ');

  _.each(params, function(param, key) {
    if (typeof param === 'object') {
      params[key] = _.map(param, function(param) {
        return parseFloat(param);
      });
    } else {
      params[key] = parseFloat(param);
    }
  });

  if (params['N'] > 3) {
    alert('N не может быть больше 3');
    return;
  }

  return params;
}

var renderChartBifurcational = function(series) {
  $('#chart').highcharts({
    chart: {
      type: 'scatter'
    },
    title: {
      text: 'Бифуркационная диаграмма'
    },
    xAxis: {
      title: {
        text: 'p'
      }
    },
    yAxis: {
      title: {
        text: 'x\'(t)'
      }
    },
    plotOptions: {
      scatter: {
        marker: {
          radius: 1
        }
      }
    },
    series: [{
      name: 'Точки диаграммы',
      data: series
    }]
  });
}

var renderChartMethod = function(seriesX, seriesXDot) {
  $('#chart').highcharts({
    chart: {
      type: 'spline'
    },
    title: {
      text: 'Исследование динамики виброударного механизма с кривошипно-шатунным возбудителем колебаний'
    },
    xAxis: {
      title: {
        text: 't'
      }
    },
    yAxis: {
      title: {
        text: 'y'
      }
    },
    plotOptions: {
      spline: {
        marker: {
          radius: 1
        }
      }
    },
    series: [{
      name: 'f(t)',
      data: seriesX
    }, {
      name: "x(t)",
      data: seriesXDot
    }]
  });
}

var getX = function(N, mu, t, t0, lambda, gamma, phi, Y0, X0, p) {
  if (N == 1) {
    return -mu * lambda[0] * (Math.cos(t) - Math.cos(t0)) -
      mu * lambda[0] * (t - t0) * Math.sin(t0) -
      p * Math.pow(t - t0, 2) / 2 +
      Y0 * (t - t0) +
      X0;
  } else if (N == 2) {
    return -mu * lambda[0] * (Math.cos(t) - Math.cos(t0)) -
      mu * (t - t0) * (lambda[0] * Math.sin(t0) + lambda[1] * gamma[0] * Math.sin(t0 - phi[0])) -
      mu * lambda[1] * gamma[0] * (Math.cos(t - phi[0]) - Math.cos(t0 - phi[0])) -
      p * ((t - t0) * (t - t0)) / 2 +
      Y0 * (t - t0) +
      X0;
  } else if (N == 3) {
    var aa = mu * (gamma[1] * lambda[1] * Math.sin(phi[1]) + gamma[2] * lambda[2] * Math.sin(phi[2]));
    var bb = mu * (lambda[0] + gamma[1] * lambda[1] * Math.cos(phi[1]) + gamma[2] * lambda[2] * Math.cos(phi[2]));
    var A = Math.sqrt(aa * aa + bb * bb);

    return -p * (t - t0) * (t - t0) / 2 +
      Y0 * (t - t0) +
      X0 -
      A * (t - t0) * Math.sin(t0 - phi[0]) -
      A * (Math.cos(t - phi[0]) - Math.cos(t0 - phi[0]));
  }
}

var getY = function(N, mu, t, t0, lambda, gamma, phi, Y0, p) {
  if (N == 1) {
    return -mu * lambda[0] * (Math.sin(t) - Math.sin(t0)) -
      p * (t - t0) +
      Y0;
  } else if (N == 2) {
    return mu * lambda[0] * (Math.sin(t) - Math.sin(t0)) +
    mu * lambda[1] * gamma[0] * (Math.sin(t - phi[0]) - Math.sin(t0 - phi[0])) -
    p * (t - t0) +
    Y0;
  } else if (N == 3) {
    var aa = mu * (gamma[1] * lambda[1] * Math.sin(phi[1]) + gamma[2] * lambda[2] * Math.sin(phi[2]));
    var bb = mu * (lambda[0] + gamma[1] * lambda[1] * Math.cos(phi[1]) + gamma[2] * lambda[2] * Math.cos(phi[2]));
    var A = Math.sqrt(aa * aa + bb * bb);
    return -p * (t - t0) +
      Y0 +
      A * (Math.sin(t - phi[0]) - Math.sin(t0 - phi[0]));
  }
}

var getFi = function(i, t, epsilon, mu, gamma, phi, gamma) {
  return epsilon[i] - mu * Math.cos(t - phi[i]) * gamma[i];
}

var getFmax = function(N, t, epsilon, mu, gamma, phi, gamma) {
  var arr = [];
  for (i = 0; i < N; i++) {
    arr.push({
      f: getFi(i, t, epsilon, mu, gamma, phi, gamma),
      i: i
    });
  }
  return _.max(arr, function(value) { return value.f; });
}

var getFpr = function(i, t, mu, gamma, phi) {
  return mu * gamma[i] * Math.sin(t - phi[i]);
}

var calculateBifurcational = function(N, R, epsilon, mu, gamma, phi, X0, Y0, t0, tn, h, lambda) {
  var series = [];

  var f_pr;
  var f_max, razn;
  var n, i;
  var t_st, t_n;
  var X, Y, arr;
  var ht;
  var t;
  var p;

  var t0_init = t0;

  for (p = 0.1; p < 2; p += h) {

    t_st = t0_init;
    ht = h;
    t = t0_init;

    // координаты неподвижной точки
    X0 = Math.PI * p * (R * (1 - 2 * lambda[0]) + 1) / ((1 + R) * (1 - lambda[0]));
    Y0 = Math.PI * p * (1 - R) / (mu * (1 + R) * (1 - lambda[0]));

    hitCounter = 0;
    calcCounter = 0;
    while (hitCounter < 100 && calcCounter < 100000) {
      calcCounter += 1;
      X = getX(N, mu, t, t0, lambda, gamma, phi, Y0, X0, p);
      Y = getY(N, mu, t, t0, lambda, gamma, phi, Y0, p);
      f_max = getFmax(N, t, epsilon, mu, gamma, phi, gamma).f;
      razn = X - f_max;
      if (razn < 0) {
        counter = 0;
        while (Math.abs(t - t_st) > 1E-3 && counter < 100) {
          counter += 1;
          ht /= 2;
          t_n = t_st + ht;
          _X = getX(N, mu, t, t0, lambda, gamma, phi, Y0, X0, p);

          f_max = getFmax(N, t, epsilon, mu, gamma, phi, gamma).f;

          if (_X - f_max < 0) {
            t = t_n;
          } else {
            t_st = t_n;
          }
        }

        f_max_temp = getFmax(N, t, epsilon, mu, gamma, phi, gamma)
        f_max = f_max_temp.f;
        f_max_i = f_max_temp.i;
        X = getX(N, mu, t, t0, lambda, gamma, phi, Y0, X0, p);

        f_pr = getFpr(f_max_i, t, mu, gamma, phi)
        Y0 = -R * Y + (1 + R) * f_pr;
        X0 = f_max;

        hitCounter += 1;
        if (hitCounter > 10) {
          series.push([2 - p, f_max]);
        }

        t0 = t;
        ht = h;
        razn_sk = Math.abs(f_pr - Y0);
        if (razn_sk <= 1E-3) {
          ht = 1E-5;
        } else {
          ht = h;
        } 
      }
      t_st = t;
      t += ht;
      t_n = 0;
      ht = h;
    }
  }

  renderChartBifurcational(series);
}

var calculate = function(N, p, R, epsilon, mu, gamma, phi, X0, Y0, t0, tn, h, lambda) {
  var f_pr;
  var f_max, razn;
  var n, i;
  var t_st, t_n;
  var X, Y, arr;

  var ht = h;

  var t_st = t0;
  var t = t0 + ht;

  var seriesX = [];
  var seriesXDot = [];

  while (t < tn) {
    X = getX(N, mu, t, t0, lambda, gamma, phi, Y0, X0, p);
    Y = getY(N, mu, t, t0, lambda, gamma, phi, Y0, p);

    f_max = getFmax(N, t, epsilon, mu, gamma, phi, gamma).f;

    seriesX.push([t, f_max]);

    razn = X - f_max;
    if (razn < 0) {

      while (Math.abs(t - t_st) > 1E-6) {
        ht /= 2;
        t_n = t_st + ht;
        _X = getX(N, mu, t, t0, lambda, gamma, phi, Y0, X0, p);

        f_max = getFmax(N, t, epsilon, mu, gamma, phi, gamma).f;

        if (_X - f_max < 0) {
          t = t_n;
        } else {
          t_st = t_n;
        }
      }

      f_max_temp = getFmax(N, t, epsilon, mu, gamma, phi, gamma)
      f_max = f_max_temp.f;
      f_max_i = f_max_temp.i;
      X = getX(N, mu, t, t0, lambda, gamma, phi, Y0, X0, p);

      f_pr = getFpr(f_max_i, t, mu, gamma, phi)
      Y0 = -R * Y + (1 + R) * f_pr;
      X0 = f_max;

      seriesXDot.push([t, X0]);

      t0 = t;
      ht = h;
      razn_sk = Math.abs(f_pr - Y0);
      if (razn_sk <= 1E-3) {
        ht = 1E-5;
      } else {
        ht = h;
      }
    } else {
      seriesXDot.push([t, X]);
    }
    t_st = t;
    t += ht;
    t_n = 0;
    ht = h;
  }

  renderChartMethod(seriesX, seriesXDot);
}