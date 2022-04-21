#include <QApplication>
#include <iostream>
#include <QLocale>
#include <QTranslator>
#include <QTextCodec>
#include <QDebug>
#include <QTime>
#include <QSplashScreen>
#include <QGLFormat>

#include "../libs/basic/logger.h"
#include "main_window.h"

#include <time.h>



class Application : public QApplication
{
public:
	Application(int& argc, char ** argv) : QApplication(argc, argv) { }
	virtual ~Application() { }

	// reimplemented from QApplication so we can throw exceptions in slots
	virtual bool notify(QObject * receiver, QEvent * event) {
		try {
			return QApplication::notify(receiver, event);
		}
		catch (std::exception& e) {
			qCritical() << "Exception thrown:" << e.what();
		}
		return false;
	}
};


int main(int argc, char **argv)
{
	srand(time(0));

#if (QT_VERSION >= QT_VERSION_CHECK(5, 6, 0))
	QApplication::setAttribute(Qt::AA_EnableHighDpiScaling);
#endif

	{
        QGLFormat format = QGLFormat::defaultFormat();
        format.setProfile(QGLFormat::CompatibilityProfile);
        format.setSampleBuffers(true); // you can also call setOption(QGL::SampleBuffers)
        format.setSamples(8);  // 8 is enough
	}

	Application app(argc, argv);

	//Locale management
	{
		//Force 'English' locale so as to get a consistent behavior everywhere
		QLocale locale = QLocale(QLocale::English);
		locale.setNumberOptions(QLocale::c().numberOptions());
		QLocale::setDefault(locale);

#ifdef Q_OS_UNIX
		//We reset the numeric locale for POSIX functions
		//See http://qt-project.org/doc/qt-5/qcoreapplication.html#locale-settings
		setlocale(LC_NUMERIC, "C");
#endif
	}

#ifndef _DEBUG
	// splash screen
	QPixmap pixmap(":/Resources/splash.png");
	QSplashScreen splash(pixmap, Qt::WindowStaysOnTopHint);
	QTime splashTimer;
	splashTimer.start();
	splash.show();
	splash.showMessage("  Starting MVStudio...");
	app.processEvents();

	//we want the splash screen to be visible a minimum amount of time (in milliseconds.)
	while (splashTimer.elapsed() < 500) {
		splash.raise();
		QApplication::processEvents(); //to let the system breath!
	}
#endif

	MainWindow window;	
	window.show();

#ifndef _DEBUG
	splash.finish(&window);
	//splash.close();
	QApplication::processEvents();
#endif

	//if (VerifyRegistry() == false)
	//	return 1;

	return app.exec();
};
