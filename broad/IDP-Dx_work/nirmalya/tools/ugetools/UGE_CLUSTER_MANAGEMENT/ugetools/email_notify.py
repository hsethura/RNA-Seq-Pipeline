import time


def email_report(message, subject, recipients, smtp_host):
    import smtplib
    from email.mime.text import MIMEText
    import getpass
    sender = getpass.getuser()
    body = message
    msg = MIMEText(body)
    msg['Subject'] = subject
    msg['To'] = ", ".join(recipients)
    msg['From'] = sender

    s = smtplib.SMTP(smtp_host)
    s.sendmail(sender, recipients, msg.as_string())
    s.quit()


class MailThrottler(object):
    def __init__(self, smtp, recipients, throttle, subject):
        self.smtp = smtp
        self.recipients = recipients
        self.throttle = throttle
        self.subject = subject
        # Set last send time in past so that first message is never throttled.
        self.last_send_time = time.time() - throttle - 1
        self.lst_messages = []

    def maybe_send(self):
        now = time.time()
        if len(self.lst_messages) > 0 and now - self.last_send_time >= self.throttle:
            email_report(message="\n\n".join(self.lst_messages), subject=self.subject,
                                      recipients=self.recipients, smtp_host=self.smtp)
            self.lst_messages = []
            self.last_send_time = now

    def add_message(self, message):
        if self.recipients is not None and len(self.recipients) > 0:
            self.lst_messages.append(message)
            self.maybe_send()